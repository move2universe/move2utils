## Tests for the `state =` parameter on mt_clean_track().

suppressPackageStartupMessages({
  library(move2); library(sf)
})

read_synthetic <- function() {
  path <- system.file("extdata", "synthetic_tracks.csv.gz",
                       package = "move2utils")
  if (nchar(path) == 0)
    path <- "inst/extdata/synthetic_tracks.csv.gz"
  d <- read.csv(gzfile(path), stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  m <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier",
    crs = 4326)
  m <- m[!sf::st_is_empty(m), ]
  m <- dplyr::arrange(m, move2::mt_track_id(m), move2::mt_time(m))
  m
}


test_that("state = NULL reproduces the no-state result", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  base <- mt_clean_track(m_A, v_max = 30, plot = FALSE,
                          remove = FALSE, silent = TRUE)
  null_state <- mt_clean_track(m_A, v_max = 30, state = NULL,
                                plot = FALSE, remove = FALSE,
                                silent = TRUE)
  expect_identical(base$is_outlier, null_state$is_outlier)
})


test_that("state vector of length nrow runs without error", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  ## Synthetic two-state assignment: split row index in half.
  st <- c(rep("rest", floor(nrow(m_A) / 2)),
          rep("flight", nrow(m_A) - floor(nrow(m_A) / 2)))
  out <- mt_clean_track(m_A, v_max = 30, state = st,
                         plot = FALSE, remove = FALSE, silent = TRUE)
  expect_s3_class(out, "move2")
  expect_equal(nrow(out), nrow(m_A))
  expect_true("is_outlier" %in% names(out))
})


test_that("state column-name dispatch works", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  m_A$beh <- c(rep("rest", floor(nrow(m_A) / 2)),
               rep("flight", nrow(m_A) - floor(nrow(m_A) / 2)))
  out <- mt_clean_track(m_A, v_max = 30, state = "beh",
                         plot = FALSE, remove = FALSE, silent = TRUE)
  expect_equal(nrow(out), nrow(m_A))
})


test_that("state errors on missing column or wrong vector length", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  expect_error(mt_clean_track(m_A, v_max = 30, state = "nope",
                              plot = FALSE, remove = FALSE, silent = TRUE),
               "not found")
  expect_error(mt_clean_track(m_A, v_max = 30, state = c("a", "b"),
                              plot = FALSE, remove = FALSE, silent = TRUE),
               "must equal nrow")
})


test_that("segments shorter than 3 fixes pass through unflagged", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  ## A 1-fix segment in the middle.
  st <- rep("flight", nrow(m_A))
  st[10] <- "rest"  # one-row segment
  ## The short segment correctly triggers the "<10 fixes" warning;
  ## that fail-safe is the behaviour we're verifying, so suppress here
  ## rather than have it surface as a test warning.
  out <- suppressWarnings(
    mt_clean_track(m_A, v_max = 30, state = st,
                   plot = FALSE, remove = FALSE, silent = TRUE))
  expect_false(out$is_outlier[10])
})


test_that("NA in state is treated as its own segment", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  st <- rep("flight", nrow(m_A))
  st[20:25] <- NA  # NA segment of 6 rows
  ## NA segment of 6 fixes is below the <10-fix guard; expected
  ## warning -- suppress so it does not surface as a test warning.
  out <- suppressWarnings(
    mt_clean_track(m_A, v_max = 30, state = st,
                   plot = FALSE, remove = FALSE, silent = TRUE))
  expect_equal(nrow(out), nrow(m_A))
})


test_that("state dispatch preserves original row order", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  ## Alternating state to maximise segment count.
  st <- rep(c("a", "b"), length.out = nrow(m_A))
  out <- mt_clean_track(m_A, v_max = 30, state = st,
                         plot = FALSE, remove = FALSE, silent = TRUE)
  ## Times should still be monotone after recombination.
  expect_true(all(diff(as.numeric(move2::mt_time(out))) >= 0))
})


test_that("state dispatch on multi-individual data works", {
  m <- read_synthetic()
  ids <- as.character(move2::mt_track_id(m))
  ## Per-fix state, alternating within each individual.
  st <- character(length(ids))
  for (id in unique(ids)) {
    i <- which(ids == id)
    st[i] <- rep_len(c("a", "b"), length(i))
  }
  out <- mt_clean_track(m, v_max = 30, state = st,
                         plot = FALSE, remove = FALSE, silent = TRUE)
  expect_equal(nrow(out), nrow(m))
})


test_that("all-stationary state subset cleans without error (degenerate prob path)", {
  ## Regression for the v0.3.3 multistate crash: a state whose subset is
  ## all-stationary (identical positions -> every step length 0) drove
  ## the probability detector's 2D histogram to a zero-height terra
  ## extent ("[rast,missing] invalid extent"), and the gap->percentile
  ## fallback on the small subset produced a zero-length flag vector.
  ## Run WITHOUT v_max so the default auto/gap prob path is exercised.
  suppressPackageStartupMessages({library(move2); library(sf)})
  set.seed(7)
  n_move <- 20L; n_rest <- 14L
  lon <- c(10 + cumsum(runif(n_move, 0.002, 0.01)), rep(10.5, n_rest))
  lat <- c(50 + cumsum(runif(n_move, 0.002, 0.01)), rep(50.5, n_rest))
  ts  <- as.POSIXct("2020-01-01", tz = "UTC") + 3600 * seq_len(n_move + n_rest)
  beh <- c(rep("move", n_move), rep("rest", n_rest))
  d <- data.frame(location.long = lon, location.lat = lat,
                  timestamp = ts, ind = "z1", beh = beh,
                  stringsAsFactors = FALSE)
  m <- move2::mt_as_move2(d, coords = c("location.long", "location.lat"),
                          time_column = "timestamp", track_id_column = "ind",
                          crs = 4326)
  out <- suppressWarnings(
    mt_clean_track(m, state = "beh", plot = FALSE,
                   remove = FALSE, silent = TRUE))
  expect_s3_class(out, "move2")
  expect_equal(nrow(out), n_move + n_rest)
  expect_true("is_outlier" %in% names(out))
  ## stationary fixes carry no kinematic signal -> must not be flagged
  expect_false(any(out$is_outlier[(n_move + 1):(n_move + n_rest)]))
})

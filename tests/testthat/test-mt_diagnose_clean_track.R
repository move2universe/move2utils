## Tests for mt_diagnose_clean_track()

suppressPackageStartupMessages({
  library(move2); library(sf)
})

read_synthetic_clean_track_result <- function() {
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
  suppressMessages(mt_clean_track(m, plot = FALSE, remove = FALSE,
                                    silent = TRUE))
}


test_that("mt_diagnose_clean_track rejects non-move2 input", {
  expect_error(mt_diagnose_clean_track(data.frame(x = 1)),
               "must be a move2 object")
})

test_that("mt_diagnose_clean_track requires the flag columns", {
  m <- move2::mt_read(system.file("extdata/synthetic_tracks.csv.gz",
                                    package = "move2utils"))
  m <- m[!sf::st_is_empty(m), ]
  expect_error(mt_diagnose_clean_track(m),
               "missing flag columns")
})

test_that("mt_diagnose_clean_track returns the expected list shape", {
  res <- read_synthetic_clean_track_result()
  diag <- suppressMessages(mt_diagnose_clean_track(res, plot = FALSE,
                                                     silent = TRUE))
  expect_type(diag, "list")
  expect_true(all(c("individual", "by_individual", "run_lengths",
                     "modes", "notes") %in% names(diag)))
  expect_s3_class(diag$by_individual, "data.frame")
  expect_true(all(c("individual", "n", "flagged", "pct") %in%
                    names(diag$by_individual)))
})

test_that("multi-individual diagnostic picks highest-flag-rate by default", {
  res <- read_synthetic_clean_track_result()
  diag <- suppressMessages(mt_diagnose_clean_track(res, plot = FALSE,
                                                     silent = TRUE))
  expect_identical(diag$individual, diag$by_individual$individual[1])
})

test_that("explicit individual = ... is honoured", {
  res <- read_synthetic_clean_track_result()
  diag <- suppressMessages(mt_diagnose_clean_track(res,
                            individual = "CPF_B",
                            plot = FALSE, silent = TRUE))
  expect_identical(diag$individual, "CPF_B")
})

test_that("unknown individual raises a clear error", {
  res <- read_synthetic_clean_track_result()
  expect_error(suppressMessages(mt_diagnose_clean_track(res,
                                  individual = "BOGUS",
                                  plot = FALSE)),
               "not found in track ids")
})

test_that("notes flag bimodality on a clearly multi-state track", {
  ## Construct a stationary-mostly track with a flight-mode burst:
  ## ~95% near-zero speeds + 5% at 5-8 m/s. Diagnostic should flag
  ## the bimodal speed distribution.
  set.seed(11)
  n <- 1000L
  is_flight <- runif(n) < 0.05
  step_m <- ifelse(is_flight, runif(n, 5, 8), abs(rnorm(n, 0.05, 0.05)))
  step_m <- pmax(step_m, 0.001)
  df <- data.frame(
    location_long = cumsum(step_m * runif(n, -1, 1)),
    location_lat  = cumsum(step_m * runif(n, -1, 1)),
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") + seq_len(n),
    individual_local_identifier = "ind",
    tag_local_identifier = "tag",
    sensor_type = "gps")
  x <- move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier", crs = 32633)
  res <- suppressMessages(
    mt_clean_track(x, plot = FALSE, remove = FALSE, silent = TRUE))
  diag <- suppressMessages(
    mt_diagnose_clean_track(res, plot = FALSE, silent = TRUE))
  expect_gte(length(diag$modes), 2L)
  expect_true(any(grepl("Panel 1.*bimodal|substantive modes",
                         diag$notes, ignore.case = TRUE)))
})

test_that(".rolling_flag_rate auto-shrinks the window on short tracks", {
  rfr <- move2utils:::.rolling_flag_rate
  ## A 3-day track with the default 7-day window previously yielded
  ## < 2 windows -> the "too few windows" placeholder (Elisa feedback
  ## 2026-06-01).  It must now produce at least `min_windows` windows.
  ts <- seq(0, 3 * 86400, by = 2 * 3600)
  io <- rep(FALSE, length(ts)); io[c(5L, 9L)] <- TRUE
  short <- rfr(ts, io, window_secs = 7 * 86400, min_windows = 8L)
  expect_gte(length(short$rate), 8L)
  expect_lt(short$window_secs, 7 * 86400)        # shrunk below requested
  ## A long track keeps the requested width untouched.
  ts2 <- seq(0, 120 * 86400, by = 3600)
  long <- rfr(ts2, rep(FALSE, length(ts2)),
              window_secs = 7 * 86400, min_windows = 8L)
  expect_equal(long$window_secs, 7 * 86400)
})

test_that("single-iteration Panel 4 renders without error (no solid-box bar)", {
  ## Converged-in-1-iteration draws a text summary, not a full-width bar.
  pf <- tempfile(fileext = ".pdf")
  grDevices::pdf(pf)
  expect_silent(
    move2utils:::.panel_iteration_cumulative(5L, "no_new_flags", NA, 1.2))
  grDevices::dev.off()
  expect_true(file.exists(pf))
})

## Helper: load a small CPF synthetic track
read_cpf <- function(track = "CPF_C") {
  d <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")),
                 stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  x <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier", crs = 4326)
  x <- x[!sf::st_is_empty(x), ]
  x <- dplyr::arrange(x, move2::mt_track_id(x), move2::mt_time(x))
  x[move2::mt_track_id(x) == track, ]
}

test_that("mt_combined_outliers rejects non-move2 input", {
  expect_error(mt_combined_outliers(NULL),
               class = "move2utils_input_not_move2")
  expect_error(mt_combined_outliers(data.frame(x = 1:3)),
               class = "move2utils_input_not_move2")
})

test_that("mt_combined_outliers attaches required columns", {
  x <- read_cpf()
  res <- suppressMessages(
    mt_combined_outliers(x, min_votes = 2, plot = FALSE))
  expect_true(all(c("is_outlier", "vote_count",
                      "flag_gap", "flag_entropy", "flag_seq") %in%
                     names(res)))
  expect_type(res$is_outlier, "logical")
  expect_type(res$vote_count, "integer")
  expect_true(all(res$vote_count %in% 0:3))
  expect_equal(nrow(res), nrow(x))
})

test_that("mt_combined_outliers min_votes gates correctly", {
  x <- read_cpf("CPF_A")
  r1 <- suppressMessages(mt_combined_outliers(x, min_votes = 1,
                                                 plot = FALSE))
  r2 <- suppressMessages(mt_combined_outliers(x, min_votes = 2,
                                                 plot = FALSE))
  r3 <- suppressMessages(mt_combined_outliers(x, min_votes = 3,
                                                 plot = FALSE))
  ## Lower min_votes => >= flags as higher min_votes
  expect_gte(sum(r1$is_outlier), sum(r2$is_outlier))
  expect_gte(sum(r2$is_outlier), sum(r3$is_outlier))
})

test_that("mt_combined_outliers plots one panel per track on multi-track input", {
  x_a <- read_cpf("CPF_A")
  x_c <- read_cpf("CPF_C")
  x_both <- rbind(x_a, x_c)

  ## Redirect the plot to a null PDF device so we can verify the call
  ## doesn't error on a multi-track input.  par() inside the function
  ## switches to mfrow = grDevices::n2mfrow(n_tracks); we restore the
  ## device's mfrow afterwards via the function's own on.exit.
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressMessages(
    mt_combined_outliers(x_both, min_votes = 2, plot = TRUE)
  ))
})

test_that("mt_combined_outliers catches known outliers on CPF_A", {
  ## CPF_A has 23 injected outliers; majority vote should catch most.
  x <- read_cpf("CPF_A")
  gt <- readRDS(system.file("extdata", "synthetic_ground_truth.rds",
                              package = "move2utils"))
  truth_idx <- gt[["CPF_A"]]$index
  truth <- rep(FALSE, nrow(x))
  truth[truth_idx] <- TRUE
  res <- suppressMessages(mt_combined_outliers(x, min_votes = 2,
                                                  plot = FALSE))
  tp <- sum(res$is_outlier & truth)
  ## catches ~14/23 (all true positives) in the canonical AEQD metric space;
  ## anchor set below that with a margin.  Was 15+ under the pre-projection
  ## flat-Earth coordinates (projection stratification 2026-06-08).
  expect_true(tp >= 13L)
})

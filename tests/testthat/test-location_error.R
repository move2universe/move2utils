## Tests for the observation-error resolver and Argos LC lookup
## (.resolve_location_error, .resolve_location_error_auto, .argos_lc_sigma).
##
## These are internal helpers used by mt_flag_outliers_bridge() to
## ingest per-fix horizontal-accuracy priors from the device.

suppressPackageStartupMessages({
  library(move2); library(sf)
})


make_track <- function(n = 30L, hacc = NULL, lc = NULL) {
  set.seed(1)
  d <- data.frame(
    id        = "T1",
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") + seq(0, by = 60, length.out = n),
    lon       = 11 + cumsum(rnorm(n, 0, 1e-3)),
    lat       = 47 + cumsum(rnorm(n, 0, 1e-3))
  )
  if (!is.null(hacc)) d$eobs_horizontal_accuracy_estimate <- hacc
  if (!is.null(lc))   d$argos_lc                          <- lc
  move2::mt_as_move2(d, coords = c("lon", "lat"),
                     time_column = "timestamp", track_id_column = "id",
                     crs = 4326)
}


# ---- .argos_lc_sigma: lookup table -----------------------------------------

test_that(".argos_lc_sigma maps standard classes to documented sigma", {
  out <- .argos_lc_sigma(c("3", "2", "1", "0", "A", "B", "Z"))
  expect_equal(out, c(250, 500, 1500, 5000, 5000, 10000, NA_real_))
})

test_that(".argos_lc_sigma is case-insensitive on letters", {
  expect_equal(.argos_lc_sigma("a"), 5000)
  expect_equal(.argos_lc_sigma("b"), 10000)
  expect_true(is.na(.argos_lc_sigma("z")))
})

test_that(".argos_lc_sigma returns NA on unknown values", {
  expect_true(is.na(.argos_lc_sigma("X")))
  expect_true(is.na(.argos_lc_sigma(NA_character_)))
  expect_true(is.na(.argos_lc_sigma("")))
})

test_that(".argos_lc_sigma handles factor input via as.character", {
  fac <- factor(c("3", "B", "1"), levels = c("3", "2", "1", "0", "A", "B", "Z"))
  expect_equal(.argos_lc_sigma(fac), c(250, 10000, 1500))
})


# ---- .resolve_location_error: NULL passes through -------------------------------

test_that(".resolve_location_error returns NULL when given NULL", {
  x <- make_track()
  expect_null(.resolve_location_error(NULL, x, nrow(x)))
})


# ---- .resolve_location_error: numeric scalar broadcasts -------------------------

test_that(".resolve_location_error broadcasts a positive scalar", {
  x <- make_track(n = 20L)
  out <- .resolve_location_error(15, x, nrow(x))
  expect_length(out, 20)
  expect_true(all(out == 15))
})

test_that(".resolve_location_error rejects negative scalar", {
  x <- make_track()
  expect_error(.resolve_location_error(-5, x, nrow(x)), "non-negative")
})

test_that(".resolve_location_error rejects NA scalar", {
  x <- make_track()
  expect_error(.resolve_location_error(NA_real_, x, nrow(x)), "non-negative")
})


# ---- .resolve_location_error: numeric vector ------------------------------------

test_that(".resolve_location_error accepts a per-fix numeric vector", {
  x <- make_track(n = 10L)
  v <- runif(10, 5, 50)
  out <- .resolve_location_error(v, x, nrow(x))
  expect_equal(out, v)
})

test_that(".resolve_location_error rejects wrong-length vector", {
  x <- make_track(n = 10L)
  expect_error(.resolve_location_error(rep(5, 7), x, nrow(x)),
                "must have length")
})

test_that(".resolve_location_error converts negative entries to NA", {
  x <- make_track(n = 5L)
  out <- .resolve_location_error(c(10, 20, -1, 30, NA), x, nrow(x))
  expect_equal(out, c(10, 20, NA_real_, 30, NA_real_))
})


# ---- .resolve_location_error: column name ---------------------------------------

test_that(".resolve_location_error reads a named column", {
  hacc <- runif(15, 5, 30)
  x <- make_track(n = 15L, hacc = hacc)
  out <- .resolve_location_error("eobs_horizontal_accuracy_estimate", x, nrow(x))
  expect_equal(out, hacc)
})

test_that(".resolve_location_error errors when column is missing", {
  x <- make_track()
  expect_error(.resolve_location_error("not_a_column", x, nrow(x)),
                "not in `x`")
})

test_that(".resolve_location_error finds dashed column variant", {
  x <- make_track(n = 8L)
  ## simulate the mt_read-on-Movebank-CSV convention with dashes
  x[["eobs-horizontal-accuracy-estimate"]] <- rep(12, 8)
  out <- .resolve_location_error("eobs_horizontal_accuracy_estimate", x, nrow(x))
  expect_equal(out, rep(12, 8))
})


# ---- .resolve_location_error: "auto" --------------------------------------------

test_that(".resolve_location_error auto picks eobs_hacc when present", {
  hacc <- runif(20, 5, 30)
  x <- make_track(n = 20L, hacc = hacc)
  out <- suppressMessages(.resolve_location_error("auto", x, nrow(x)))
  expect_equal(out, hacc)
})

test_that(".resolve_location_error auto falls back to argos_lc", {
  lc <- c("3", "2", "1", "B", "A", "0", "3", "2", "Z", "B")
  x  <- make_track(n = 10L, lc = lc)
  out <- suppressMessages(.resolve_location_error("auto", x, nrow(x)))
  expected <- c(250, 500, 1500, 10000, 5000, 5000, 250, 500, NA_real_, 10000)
  expect_equal(out, expected)
})

test_that(".resolve_location_error auto returns NULL when neither column present", {
  x <- make_track()
  out <- suppressMessages(.resolve_location_error("auto", x, nrow(x)))
  expect_null(out)
})


# ---- .resolve_location_error: invalid input -------------------------------------

test_that(".resolve_location_error rejects bogus types", {
  x <- make_track()
  expect_error(.resolve_location_error(list(a = 1), x, nrow(x)),
                "must be NULL")
  expect_error(.resolve_location_error(c("a", "b"), x, nrow(x)),
                "must be NULL")
})

## Tests for mt_suggest_dbbmm_window().

suppressPackageStartupMessages({
  library(move2); library(sf)
})


make_regular <- function(n = 50L, dt_secs = 600) {
  d <- data.frame(
    id        = "T1",
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") + (seq_len(n) - 1L) * dt_secs,
    lon       = 11 + cumsum(rnorm(n, 0, 1e-3)),
    lat       = 47 + cumsum(rnorm(n, 0, 1e-3))
  )
  move2::mt_as_move2(d, coords = c("lon", "lat"),
                     time_column = "timestamp", track_id_column = "id",
                     crs = 4326)
}


# ---- input validation ------------------------------------------------------

test_that("mt_suggest_dbbmm_window rejects non-move2 input", {
  expect_error(mt_suggest_dbbmm_window(data.frame(x = 1)),
               class = "move2utils_input_not_move2")
})

test_that("mt_suggest_dbbmm_window rejects bad target_hours", {
  m <- make_regular()
  expect_error(mt_suggest_dbbmm_window(m, target_hours = 0),
               class = "move2utils_mt_suggest_dbbmm_window_bad_target_hours")
  expect_error(mt_suggest_dbbmm_window(m, target_hours = -1),
               class = "move2utils_mt_suggest_dbbmm_window_bad_target_hours")
  expect_error(mt_suggest_dbbmm_window(m, target_hours = c(1, 2)),
               class = "move2utils_mt_suggest_dbbmm_window_bad_target_hours")
})


# ---- structure of return value --------------------------------------------

test_that("mt_suggest_dbbmm_window returns expected fields", {
  m <- make_regular(n = 100L, dt_secs = 600)
  s <- suppressMessages(mt_suggest_dbbmm_window(m, target_hours = 4, plot = FALSE))
  expect_s3_class(s, "mt_dbbmm_window_suggestion")
  for (nm in c("window_size", "margin", "median_dt_secs",
               "target_hours", "window_hours", "n", "n_windows")) {
    expect_true(nm %in% names(s),
                info = sprintf("missing field: %s", nm))
  }
  ## window_size and margin must be odd integers
  expect_true(s$window_size %% 2L == 1L)
  expect_true(s$margin      %% 2L == 1L)
  expect_true(s$window_size >= 11L)
  expect_true(s$margin      >= 3L)
  expect_true(2L * s$margin + 1L <= s$window_size)
})


# ---- target_hours scaling --------------------------------------------------

test_that("larger target_hours gives larger window_size", {
  m <- make_regular(n = 500L, dt_secs = 600)
  s_small <- suppressMessages(
    mt_suggest_dbbmm_window(m, target_hours = 1, plot = FALSE))
  s_big   <- suppressMessages(
    mt_suggest_dbbmm_window(m, target_hours = 8, plot = FALSE))
  expect_lt(s_small$window_size, s_big$window_size)
})

test_that("window_size is bounded above by floor(n/4)", {
  m <- make_regular(n = 60L, dt_secs = 600)  # 60 fixes, n/4 = 15
  s <- suppressMessages(
    mt_suggest_dbbmm_window(m, target_hours = 100, plot = FALSE))
  expect_lte(s$window_size, 15L)
})


# ---- floor at window_size = 11 --------------------------------------------

test_that("very small target_hours hits the window_size = 11 floor", {
  m <- make_regular(n = 50L, dt_secs = 60)  # 1 min sampling
  s <- suppressMessages(
    mt_suggest_dbbmm_window(m, target_hours = 0.05, plot = FALSE))
  expect_equal(s$window_size, 11L)
})


# ---- short track returns NULL ----------------------------------------------

test_that("track with fewer than 11 fixes returns NULL", {
  m <- make_regular(n = 8L)
  s <- suppressMessages(mt_suggest_dbbmm_window(m, plot = FALSE))
  expect_null(s)
})


# ---- multi-track dispatch --------------------------------------------------

test_that("multi-track input returns one suggestion per individual", {
  d <- rbind(
    data.frame(
      id = "T1",
      timestamp = as.POSIXct("2024-01-01", tz = "UTC") +
        (seq_len(100L) - 1L) * 600,
      lon = 11 + cumsum(rnorm(100L, 0, 1e-3)),
      lat = 47 + cumsum(rnorm(100L, 0, 1e-3))),
    data.frame(
      id = "T2",
      timestamp = as.POSIXct("2024-02-01", tz = "UTC") +
        (seq_len(50L) - 1L) * 1800,
      lon = 12 + cumsum(rnorm(50L, 0, 1e-3)),
      lat = 48 + cumsum(rnorm(50L, 0, 1e-3))))
  m <- move2::mt_as_move2(d, coords = c("lon", "lat"),
                            time_column = "timestamp",
                            track_id_column = "id", crs = 4326)
  out <- suppressMessages(mt_suggest_dbbmm_window(m, target_hours = 4, plot = FALSE))
  expect_named(out, c("T1", "T2"))
  expect_s3_class(out$T1, "mt_dbbmm_window_suggestion")
  expect_s3_class(out$T2, "mt_dbbmm_window_suggestion")
})


# ---- print method ---------------------------------------------------------

test_that("print method prints expected fields", {
  m <- make_regular(n = 100L, dt_secs = 600)
  s <- suppressMessages(mt_suggest_dbbmm_window(m, plot = FALSE))
  out <- capture.output(print(s))
  expect_true(any(grepl("window_size", out)))
  expect_true(any(grepl("margin",      out)))
  expect_true(any(grepl("n locations", out)))
})

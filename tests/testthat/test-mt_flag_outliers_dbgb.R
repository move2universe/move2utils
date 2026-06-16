## Tests for mt_flag_outliers_dbgb().
##
## Coverage focus: input validation, structural output (columns +
## attributes), edge cases (short tracks), multi-track dispatch,
## auto-projection of lon/lat input, and the dBGB variance branch.
## Behavioural detection on curated outlier data is left to the
## integration tests in tests/testthat/test-mt_clean_track.R, which
## exercises the same primitive at the cleaning-pipeline level.

skip_if_not_installed("move2")
skip_if_not_installed("sf")


## ---- helpers --------------------------------------------------------

.make_leroy_proj <- function(n = 80) {
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  leroy   <- fishers[move2::mt_track_id(fishers) == "M4", ][seq_len(n), ]
  sf::st_transform(leroy, move2::mt_aeqd_crs(leroy))
}

.make_two_track_proj <- function(n = 60) {
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  ids     <- unique(move2::mt_track_id(fishers))[1:2]
  parts   <- lapply(ids, function(id)
    fishers[move2::mt_track_id(fishers) == id, ][seq_len(n), ])
  sf::st_transform(do.call(rbind, parts),
                   move2::mt_aeqd_crs(parts[[1]]))
}


## ---- tests ----------------------------------------------------------

test_that("rejects non-move2 input", {
  expect_error(mt_flag_outliers_dbgb(data.frame(x = 1:5, y = 1:5)),
               "must be a move2 object")
})

test_that("rejects malformed pre_peel_v_max", {
  leroy <- .make_leroy_proj(60)
  expect_error(suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(leroy, pre_peel_v_max = -1,
                           plot = FALSE, silent = TRUE))),
    "positive scalar")
  expect_error(suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(leroy, pre_peel_v_max = c(5, 10),
                           plot = FALSE, silent = TRUE))),
    "positive scalar")
  expect_error(suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(leroy, pre_peel_v_max = NA_real_,
                           plot = FALSE, silent = TRUE))),
    "positive scalar")
})

test_that("rejects malformed z_threshold and residual_max", {
  leroy <- .make_leroy_proj(60)
  ## Numeric z_threshold is only meaningful with a parametric
  ## threshold method; the gap method is data-driven.
  expect_error(suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(leroy, z_threshold = -1,
                           z_threshold_method = "bonferroni",
                           plot = FALSE, silent = TRUE))),
    "positive scalar")
  expect_error(suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(leroy, residual_max = "huge",
                           plot = FALSE, silent = TRUE))),
    "positive scalar, Inf, or NULL")
})

test_that("returns expected columns and attributes on a single track", {
  ## Default is now variance = "dbgb" + z_threshold_method = "gap"
  ## (envelope rule: per-channel gap thresholds on Z_para, Z_orth,
  ## Z_chisq, with bridge_z_class labelling which channel(s) fired).
  leroy <- .make_leroy_proj(80)
  out   <- suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(leroy, plot = FALSE, silent = TRUE)))

  expected_cols <- c("bridge_residual", "bridge_width",
                     "sigma2_motion_para", "sigma2_motion_orth",
                     "bridge_z", "bridge_z_chisq", "bridge_z_class",
                     "bridge_z_para", "bridge_z_orth",
                     "bridge_iteration", "is_outlier")
  for (col in expected_cols) expect_true(col %in% names(out),
                                          info = paste("missing:", col))

  expect_type(out$is_outlier, "logical")
  expect_equal(length(out$is_outlier), nrow(out))
  expect_equal(nrow(out), nrow(leroy))

  ## attributes set by the function
  expect_true(!is.null(attr(out, "z_thresholds")))
  expect_true(!is.null(attr(out, "z_threshold")))
  expect_true(!is.null(attr(out, "z_threshold_method")))
  expect_true(!is.null(attr(out, "residual_max")))
  expect_true(!is.null(attr(out, "window_size")))
  expect_true(!is.null(attr(out, "location_error")))
})

test_that("variance = 'dbbmm' (legacy isotropic) returns scalar-Z schema", {
  leroy <- .make_leroy_proj(80)
  out   <- suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(leroy, variance = "dbbmm",
                           plot = FALSE, silent = TRUE)))
  expected_cols <- c("bridge_residual", "bridge_width", "sigma2_motion",
                     "bridge_z", "bridge_z_chisq", "bridge_z_class",
                     "bridge_iteration", "is_outlier")
  for (col in expected_cols) expect_true(col %in% names(out),
                                          info = paste("missing:", col))
  ## per-axis channels are dBGB-only; absent under dBBMM
  expect_false("bridge_z_para" %in% names(out))
  expect_false("bridge_z_orth" %in% names(out))
})

test_that("short track (< 10 locations) returns NA columns and FALSE flags", {
  leroy <- .make_leroy_proj(80)[1:8, ]
  out   <- suppressWarnings(suppressMessages(
    mt_flag_outliers_dbgb(leroy, plot = FALSE, silent = TRUE)))

  expect_true(all(is.na(out$bridge_residual)))
  expect_true(all(is.na(out$bridge_z)))
  expect_true(all(out$is_outlier == FALSE))
})

test_that("multi-track input dispatches per individual and rebinds rows", {
  small <- .make_two_track_proj(60)
  out   <- suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(small, plot = FALSE, silent = TRUE)))

  expect_equal(nrow(out), nrow(small))
  expect_true("is_outlier" %in% names(out))
  ## both individuals should be represented in the output
  expect_setequal(unique(move2::mt_track_id(out)),
                  unique(move2::mt_track_id(small)))
})

test_that("lon/lat input is auto-projected and round-tripped", {
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  leroy   <- fishers[move2::mt_track_id(fishers) == "M4", ][seq_len(80), ]
  expect_true(sf::st_is_longlat(leroy))

  out <- suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(leroy, plot = FALSE, silent = TRUE)))

  expect_true(sf::st_is_longlat(out))
  expect_true("is_outlier" %in% names(out))
  expect_equal(nrow(out), nrow(leroy))
})

test_that("z_threshold_method = 'bonferroni' is the parametric envelope (per-channel)", {
  ## Parametric mode applies per-channel Bonferroni-Z thresholds with
  ## joint FWER <= alpha (alpha/2 chisq + alpha/4 each per-axis).
  ## All three channels emit finite thresholds under variance="dbgb".
  leroy <- .make_leroy_proj(80)
  out   <- suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(leroy, variance = "dbgb",
                           z_threshold_method = "bonferroni",
                           plot = FALSE, silent = TRUE)))
  thr <- attr(out, "z_thresholds")
  expect_true(is.list(thr))
  expect_true(is.numeric(thr$para)  && is.finite(thr$para)  && thr$para  > 0)
  expect_true(is.numeric(thr$orth)  && is.finite(thr$orth)  && thr$orth  > 0)
  expect_true(is.numeric(thr$chisq) && is.finite(thr$chisq) && thr$chisq > 0)
})

test_that("variance = 'dbbmm' parametric mode keeps per-axis thresholds NA", {
  ## Under dBBMM there is no per-axis decomposition; only the chisq
  ## channel exists.
  leroy <- .make_leroy_proj(80)
  out   <- suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(leroy, variance = "dbbmm",
                           z_threshold_method = "bonferroni",
                           plot = FALSE, silent = TRUE)))
  thr <- attr(out, "z_thresholds")
  expect_true(is.na(thr$para))
  expect_true(is.na(thr$orth))
  expect_true(is.finite(thr$chisq) && thr$chisq > 0)
})

test_that("remove = TRUE returns only non-flagged rows", {
  leroy <- .make_leroy_proj(80)
  full  <- suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(leroy, plot = FALSE, silent = TRUE,
                           remove = FALSE)))
  trimmed <- suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(leroy, plot = FALSE, silent = TRUE,
                           remove = TRUE)))

  expect_equal(nrow(trimmed), sum(!full$is_outlier))
})

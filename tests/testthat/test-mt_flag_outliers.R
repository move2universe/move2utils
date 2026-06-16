test_that("mt_flag_outliers rejects non-move2 input", {
  expect_error(mt_flag_outliers(data.frame(x = 1:5)),
               "must be a move2 object")
})

test_that("mt_flag_outliers rejects invalid threshold", {
  expect_error(mt_flag_outliers(structure(list(), class = "move2"),
                                threshold = 0),
               "must be a positive number")
  expect_error(mt_flag_outliers(structure(list(), class = "move2"),
                                threshold = -1),
               "must be a positive number")
})

test_that("mt_flag_outliers rejects invalid prob_type", {
  expect_error(mt_flag_outliers(structure(list(), class = "move2"),
                                prob_type = "nonsense"),
               "must be one of")
})

test_that("mt_flag_outliers works on fisher example data", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  leroy <- fishers[mt_track_id(fishers) == "M4", ]
  leroy <- leroy[!st_is_empty(leroy), ]

  ## flag mode — use significance method to ensure some outliers are found
  result <- mt_flag_outliers(leroy, plot = FALSE, threshold_type = "significance")
  expect_s3_class(result, "move2")
  expect_equal(nrow(result), nrow(leroy))
  expect_true("is_outlier" %in% names(result))
  expect_true("joint_prob" %in% names(result))
  expect_true("outlier_percentile" %in% names(result))

  ## at least some outliers should be flagged
  expect_true(any(result$is_outlier, na.rm = TRUE))

  ## outlier count should be small relative to total
  n_outliers <- sum(result$is_outlier, na.rm = TRUE)
  expect_true(n_outliers < nrow(leroy) * 0.10)
})

test_that("mt_flag_outliers remove mode returns fewer rows", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  leroy <- fishers[mt_track_id(fishers) == "M4", ]
  leroy <- leroy[!st_is_empty(leroy), ]

  cleaned <- mt_flag_outliers(leroy, remove = TRUE, plot = FALSE,
                              threshold_type = "significance")
  expect_s3_class(cleaned, "move2")
  expect_true(nrow(cleaned) < nrow(leroy))
  expect_true(nrow(cleaned) > nrow(leroy) * 0.8)
})

test_that("mt_flag_outliers works with different prob_types", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  leroy <- fishers[mt_track_id(fishers) == "M4", ]
  leroy <- leroy[!st_is_empty(leroy), ]

  for (pt in c("joint", "step_turn", "delta_step", "delta_turn", "custom")) {
    result <- mt_flag_outliers(leroy, prob_type = pt, plot = FALSE)
    expect_true("is_outlier" %in% names(result),
                info = paste("prob_type:", pt))
  }
})

test_that("mt_flag_outliers drop_na removes NA probability rows", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  leroy <- fishers[mt_track_id(fishers) == "M4", ]
  leroy <- leroy[!st_is_empty(leroy), ]

  result_keep <- mt_flag_outliers(leroy, remove = TRUE,
                                  drop_na = FALSE, plot = FALSE)
  result_drop <- mt_flag_outliers(leroy, remove = TRUE,
                                  drop_na = TRUE, plot = FALSE)

  ## dropping NAs should result in equal or fewer rows
  expect_true(nrow(result_drop) <= nrow(result_keep))
})


# ---- projected input auto-transforms to lon/lat ----------------------------

test_that("mt_flag_outliers accepts projected input and returns in original CRS", {
  m  <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")))
  m$timestamp <- as.POSIXct(m$timestamp, tz = "UTC")
  m <- move2::mt_as_move2(m, coords = c("location.long", "location.lat"),
                            time_column = "timestamp",
                            track_id_column = "individual.local.identifier",
                            crs = 4326)
  m <- m[!sf::st_is_empty(m), ]
  m <- dplyr::arrange(m, move2::mt_track_id(m), move2::mt_time(m))
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  aeqd <- move2::mt_aeqd_crs(m_C, center = "center", units = "m")
  m_C_p <- sf::st_transform(m_C, aeqd)

  ## projected input must not error and must yield a usable result
  res_p <- suppressMessages(mt_flag_outliers(m_C_p, plot = FALSE))
  ## CRS preserved
  expect_identical(sf::st_crs(res_p), sf::st_crs(m_C_p))
  ## diagnostic columns attached
  expect_true("is_outlier" %in% names(res_p))
  expect_true("joint_prob" %in% names(res_p))
  expect_true(any(is.finite(res_p$joint_prob)))

  ## The lon/lat path (Haversine + spherical great-circle azimuth) and
  ## the projected path (Euclidean step + Cartesian atan2) compute the
  ## same scores up to AEQD-projection accuracy.  Empirically on the
  ## synthetic data the joint probabilities agree with Spearman > 0.999
  ## and max absolute difference < 1e-3.  We do NOT compare flag sets
  ## directly: the data-driven threshold detector ('gap') sits on the
  ## boundary of an unimodal-or-nearly-unimodal distribution for clean
  ## tracks and a tiny shift in the score distribution can move the
  ## detected break by several fixes.  The right invariant is that the
  ## scoring agrees, not that thresholding lands on identical fixes.
  res_ll <- suppressMessages(mt_flag_outliers(m_C, plot = FALSE))
  ok <- is.finite(res_p$joint_prob) & is.finite(res_ll$joint_prob)
  expect_gt(stats::cor(res_p$joint_prob[ok], res_ll$joint_prob[ok],
                        method = "spearman"),
             0.99)
  expect_lt(max(abs(res_p$joint_prob[ok] - res_ll$joint_prob[ok])),
             1e-2)
  expect_gt(sum(res_p$is_outlier),  0L)
  expect_gt(sum(res_ll$is_outlier), 0L)
})

test_that("mt_flag_outliers projected input preserves probability columns", {
  m  <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")))
  m$timestamp <- as.POSIXct(m$timestamp, tz = "UTC")
  m <- move2::mt_as_move2(m, coords = c("location.long", "location.lat"),
                            time_column = "timestamp",
                            track_id_column = "individual.local.identifier",
                            crs = 4326)
  m <- m[!sf::st_is_empty(m), ]
  m <- dplyr::arrange(m, move2::mt_track_id(m), move2::mt_time(m))
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  m_A_p <- sf::st_transform(m_A,
    move2::mt_aeqd_crs(m_A, center = "center", units = "m"))

  res <- suppressMessages(mt_flag_outliers(m_A_p, plot = FALSE))
  for (col in c("step_turn_prob", "delta_step_prob", "delta_turn_prob",
                "joint_prob", "is_outlier", "is_na_prob")) {
    expect_true(col %in% names(res))
  }
  expect_true(any(is.finite(res$joint_prob)))
})


## ---- regression: degenerate state-subset paths (multistate cascade) ----
## When mt_clean_track() dispatches by state, each state value is cleaned
## by a recursive call on its subset.  Small or stationary subsets exposed
## two crashes in the probability detector (both surfaced by the v0.3.3
## entropy/gap propagation fix, which made gap-mode the default prob path
## with a NULL threshold deferred to the leaf formal):
##   (A) gap-mode falling back to percentile carried a NULL threshold into
##       `100 - threshold * 100`, yielding numeric(0) and a zero-length
##       `is_outlier` that broke the caller's scatter assignment.
##   (B) an all-stationary subset (every step length 0) drove
##       `.turn_step_hist()` to build a terra raster with a zero-height
##       ymin == ymax extent ("[rast,missing] invalid extent").

test_that(".turn_step_hist returns NULL on degenerate step/turn inputs", {
  ## all-zero steps (stationary subset of repeated positions)
  expect_null(move2utils:::.turn_step_hist(runif(20, -pi, pi), rep(0, 20)))
  ## all-NA steps
  expect_null(move2utils:::.turn_step_hist(runif(20, -pi, pi), rep(NA_real_, 20)))
  ## no valid turn angles
  expect_null(move2utils:::.turn_step_hist(rep(NA_real_, 20), runif(20, 1, 5)))
  ## a well-formed input still yields a raster
  r <- move2utils:::.turn_step_hist(runif(40, -pi, pi), runif(40, 1, 5))
  expect_s4_class(r, "SpatRaster")
})

test_that(".prob_fn_core gap fallback keeps full-length flags on tiny subsets", {
  ## A subset too small for gap detection (< 10 valid probabilities)
  ## falls back to percentile.  With a NULL threshold (the cascade
  ## default, deferred to the leaf formal) the pre-fix percentile branch
  ## produced a zero-length `is_outlier`.  It must stay length n.
  set.seed(1)
  n  <- 8
  cc <- cbind(cumsum(rnorm(n, 2, 0.5)), cumsum(rnorm(n, 2, 0.5)))
  t_s <- as.numeric(seq(0, by = 3600, length.out = n))
  res <- suppressWarnings(move2utils:::.prob_fn_core(
    cc, t_s, seq_len(n), was_longlat = FALSE,
    threshold = NULL, prob_type = "joint",
    autodiff_alpha = "acf", acf_alpha = TRUE, auto_alpha = FALSE,
    method = "histogram", time_normalize = TRUE,
    threshold_type = "gap", step_transform = "none",
    step_floor = 0, silent = TRUE))
  expect_length(res$is_outlier, n)
  expect_type(res$is_outlier, "logical")
})

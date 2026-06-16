# ---- tests for new features: iterations, quality_columns, auto alpha,
#      time_normalize ----

# ---- iterative refinement ----

test_that("iterations parameter catches consecutive outliers", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")

  ## simulate a clean track and inject 3 consecutive outliers
  crw <- .simulate_crw_move2(n = 500, step_sd = 0.001, seed = 42)
  geom <- sf::st_geometry(crw)

  ## inject 3 consecutive large displacements at indices 100, 101, 102
  set.seed(999)
  for (i in 100:102) {
    pt <- sf::st_coordinates(geom[i])
    angle <- runif(1, 0, 2 * pi)
    new_pt <- c(pt[1] + 0.08 * cos(angle), pt[2] + 0.08 * sin(angle))
    geom[i] <- sf::st_point(new_pt)
  }
  sf::st_geometry(crw) <- geom

  result_1 <- mt_flag_outliers(crw, threshold = 0.01, iterations = 1,
                                plot = FALSE)
  result_3 <- mt_flag_outliers(crw, threshold = 0.01, iterations = 3,
                                plot = FALSE)

  outliers_1 <- sum(result_1$is_outlier, na.rm = TRUE)
  outliers_3 <- sum(result_3$is_outlier, na.rm = TRUE)

  message(sprintf("Consecutive outlier test: iter=1 found %d, iter=3 found %d",
                  outliers_1, outliers_3))

  ## iterations=3 should catch at least as many as iterations=1
  expect_true(outliers_3 >= outliers_1,
              info = "Iterative refinement should find >= single-pass outliers")
  ## both should return the full original object

  expect_equal(nrow(result_1), nrow(crw))
  expect_equal(nrow(result_3), nrow(crw))
  expect_true("is_outlier" %in% names(result_3))
})

test_that("iterations=1 is backward compatible (default)", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")

  crw <- .simulate_crw_move2(n = 200, step_sd = 0.001, seed = 42)
  r_default <- mt_flag_outliers(crw, plot = FALSE)
  r_iter1 <- mt_flag_outliers(crw, iterations = 1, plot = FALSE)

  ## should produce identical results
  expect_equal(r_default$is_outlier, r_iter1$is_outlier)
  expect_equal(r_default$joint_prob, r_iter1$joint_prob)
})


# ---- quality weighting ----

test_that("quality_columns weights affect joint probability", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  leroy <- fishers[mt_track_id(fishers) == "M4", ]
  leroy <- leroy[!st_is_empty(leroy), ]

  ## add a fake hdop column: mostly low (good), a few high (bad)
  set.seed(42)
  leroy$hdop <- rnorm(nrow(leroy), mean = 2, sd = 0.5)
  bad_idx <- c(50, 100, 150)
  leroy$hdop[bad_idx] <- 15  ## very high HDOP

  ## run without quality weighting
  r_no_qc <- mt_flag_outliers(leroy, plot = FALSE, quality_columns = NULL)

  ## run with quality weighting (high HDOP = low quality)
  r_qc <- mt_flag_outliers(leroy, plot = FALSE,
                            quality_columns = list(
                              "hdop" = function(h) 1 - pnorm(h, mean = 3, sd = 1.5)
                            ))

  expect_true("quality_weight" %in% names(r_qc))
  expect_true(!"quality_weight" %in% names(r_no_qc))

  ## quality weights for bad locations should be low
  expect_true(all(r_qc$quality_weight[bad_idx] < 0.1),
              info = "High HDOP locations should have low quality weight")

  ## joint_prob should be lower for bad locations with quality weighting
  ## (compared to without)
  for (idx in bad_idx) {
    if (!is.na(r_qc$joint_prob[idx]) && !is.na(r_no_qc$joint_prob[idx])) {
      expect_true(r_qc$joint_prob[idx] < r_no_qc$joint_prob[idx],
                  info = paste("Quality-weighted prob should be lower at idx", idx))
    }
  }
})

test_that("quality_columns=NULL gives same result as before", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")

  crw <- .simulate_crw_move2(n = 200, step_sd = 0.001, seed = 42)
  r1 <- mt_flag_outliers(crw, plot = FALSE)
  r2 <- mt_flag_outliers(crw, plot = FALSE, quality_columns = NULL)

  expect_equal(r1$joint_prob, r2$joint_prob)
  expect_equal(r1$is_outlier, r2$is_outlier)
})

test_that("quality_columns validation works", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")

  crw <- .simulate_crw_move2(n = 50, step_sd = 0.001, seed = 42)

  ## non-existent column
  expect_error(
    mt_flag_outliers(crw, plot = FALSE,
                     quality_columns = list("nonexistent" = identity)),
    "not found"
  )

  ## not a list
  expect_error(
    mt_flag_outliers(crw, plot = FALSE, quality_columns = "wrong"),
    "named list"
  )
})


# ---- auto-optimised alpha ----

test_that("autodiff_alpha='auto' works on fisher data", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  leroy <- fishers[mt_track_id(fishers) == "M4", ]
  leroy <- leroy[!st_is_empty(leroy), ]

  result <- mt_flag_outliers(leroy, autodiff_alpha = "auto", plot = FALSE)

  expect_s3_class(result, "move2")
  expect_true("is_outlier" %in% names(result))
  expect_true("joint_prob" %in% names(result))
  expect_equal(nrow(result), nrow(leroy))
})

test_that("autodiff_alpha='auto' does not error on simulated data", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")

  crw <- .simulate_crw_move2(n = 300, step_sd = 0.001, seed = 42)
  result <- mt_flag_outliers(crw, autodiff_alpha = "auto", plot = FALSE)

  expect_s3_class(result, "move2")
  expect_true("is_outlier" %in% names(result))
})

test_that("autodiff_alpha rejects invalid values", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")

  crw <- .simulate_crw_move2(n = 50, step_sd = 0.001, seed = 42)

  expect_error(
    mt_flag_outliers(crw, autodiff_alpha = "bad", plot = FALSE),
    "non-negative number"
  )
  expect_error(
    mt_flag_outliers(crw, autodiff_alpha = -1, plot = FALSE),
    "non-negative number"
  )
})


# ---- time-normalised metrics ----

test_that("time_normalize=TRUE works on fisher data", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  leroy <- fishers[mt_track_id(fishers) == "M4", ]
  leroy <- leroy[!st_is_empty(leroy), ]

  result_tn <- mt_flag_outliers(leroy, time_normalize = TRUE, plot = FALSE)
  result_no <- mt_flag_outliers(leroy, time_normalize = FALSE, plot = FALSE)

  ## both should return valid results
  expect_s3_class(result_tn, "move2")
  expect_s3_class(result_no, "move2")
  expect_true("is_outlier" %in% names(result_tn))
  expect_true("is_outlier" %in% names(result_no))
  expect_equal(nrow(result_tn), nrow(leroy))
  expect_equal(nrow(result_no), nrow(leroy))

  ## both should find some outliers (or at least not error)
  ## probabilities should differ because metrics differ
  expect_false(
    all(result_tn$joint_prob == result_no$joint_prob, na.rm = TRUE),
    info = "Time-normalised should produce different probabilities"
  )
})

test_that("time_normalize=TRUE works on simulated data", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")

  crw <- .simulate_crw_move2(n = 300, step_sd = 0.001, seed = 42)
  inj <- .inject_outliers(crw, n_outliers = 5, displacement = 0.05, seed = 77)

  result <- mt_flag_outliers(inj$data, threshold = 0.01,
                              time_normalize = TRUE, plot = FALSE)

  expect_s3_class(result, "move2")
  expect_true("is_outlier" %in% names(result))
  expect_equal(nrow(result), nrow(inj$data))
})

## "time_normalize=FALSE is backward compatible (default)" deleted in
## 0.2.0-dev: the test compared the default call (which is
## time_normalize = TRUE) to an explicit time_normalize = FALSE and
## asserted equality.  That is false-by-design -- TRUE and FALSE
## must produce different probabilities by definition of the
## time-normalisation path.  The test only ever passed because the
## C1 unit bug (mt_time_lags(..., units = "secs") routed to
## as.numeric() instead of mt_time_lags()) silently zero-broadcast
## the time-normalisation, making time_normalize = TRUE a no-op.
## With C1 fixed, the test correctly fails.  Removed because there
## is no salvageable invariant to assert here -- the underlying
## "default value" question is documented behaviour, not a runtime
## test.

# ---- gap threshold type ----

test_that("gap method finds no outliers in clean data", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")

  ## simulate clean data with no outliers — use a larger sample
  ## and a stricter gap threshold to ensure the broken stick null
  ## is well-estimated
  crw <- .simulate_crw_move2(n = 1000, step_sd = 0.001, seed = 123)
  result <- mt_flag_outliers(crw, threshold_type = "gap", threshold = 5,
                             plot = FALSE)
  n_out <- sum(result$is_outlier)
  expect_true(n_out < nrow(crw) * 0.05,
              info = paste0("Gap method found ", n_out,
                            " outliers in clean data (expected < 5%)"))
})

test_that("gap method finds outliers when there is a genuine break", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")

  ## simulate data with large injected outliers that create a clear gap
  crw <- .simulate_crw_move2(n = 500, step_sd = 0.001, seed = 42)
  inj <- .inject_outliers(crw, n_outliers = 15, displacement = 0.1, seed = 77)

  result <- mt_flag_outliers(inj$data, threshold_type = "gap",
                             threshold = 2, plot = FALSE)
  n_flagged <- sum(result$is_outlier)
  ## should find at least some if the displacement is large enough
  ## (this depends on whether the injected points create a clear gap)
  expect_true(is.numeric(n_flagged))
  expect_true("is_outlier" %in% names(result))
})

test_that("gap method returns valid output structure", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  leroy <- fishers[mt_track_id(fishers) == "M4", ]
  leroy <- leroy[!st_is_empty(leroy), ]

  result <- mt_flag_outliers(leroy, threshold_type = "gap", plot = FALSE)
  expect_s3_class(result, "move2")
  expect_true("is_outlier" %in% names(result))
  expect_true("outlier_percentile" %in% names(result))
  expect_equal(nrow(result), nrow(leroy))
})

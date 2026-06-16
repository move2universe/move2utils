## Regression tests for the move->move2 porting bug whereby a per-row
## location_error vector passed to a multi-track move2 was not sliced
## per-track before being handed to the single-track variance helper.
## See NEWS entry for v0.3.1 dev cycle.

.make_fisher_multi <- function(min_n = 50) {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]

  ## Project to AEQD so dBBMM accepts the object.
  fishers_proj <- st_transform(fishers, mt_aeqd_crs(fishers))

  ## Keep only tracks long enough for the default window.
  ok_ids <- names(which(table(mt_track_id(fishers_proj)) >= min_n))
  fishers_proj[mt_track_id(fishers_proj) %in% ok_ids, ]
}

test_that("mt_dbbmm_variance accepts a per-fix vector on multi-track input", {
  x <- .make_fisher_multi()
  loc_err <- rep(25, nrow(x))

  res <- suppressMessages(suppressWarnings(
    mt_dbbmm_variance(x, location_error = loc_err,
                       window_size = 31, margin = 11)
  ))
  expect_type(res, "list")
  expect_true(all(vapply(res, inherits, logical(1), "mt_dbbmm_variance")))
  expect_setequal(names(res), as.character(unique(mt_track_id(x))))
})

test_that("multi-track vector path is byte-identical to per-track scalar fit", {
  x <- .make_fisher_multi()
  ids <- as.character(mt_track_id(x))

  ## Scalar path (existing, byte-identical baseline).
  scalar_list <- mt_dbbmm_variance(x, location_error = 25,
                                    window_size = 31, margin = 11)

  ## Per-row vector path (the post-fix regression).
  vec_list <- mt_dbbmm_variance(x, location_error = rep(25, nrow(x)),
                                 window_size = 31, margin = 11)

  expect_setequal(names(scalar_list), names(vec_list))
  for (nm in names(scalar_list)) {
    expect_equal(scalar_list[[nm]]$variance,
                 vec_list[[nm]]$variance,
                 tolerance = .Machine$double.eps^0.5,
                 info = paste("track", nm))
  }
})

test_that("per-track slice from multi-track call matches a direct single-track fit", {
  x <- .make_fisher_multi()
  ids <- as.character(mt_track_id(x))

  ## Inject a track-specific error so we can detect mis-slicing: each
  ## track gets a different uniform sigma.
  uid <- unique(ids)
  sigma_by_id <- setNames(seq(10, 10 + 5 * (length(uid) - 1), by = 5), uid)
  loc_err <- sigma_by_id[ids]

  multi <- mt_dbbmm_variance(x, location_error = loc_err,
                              window_size = 31, margin = 11)

  for (id in names(multi)) {
    sub <- x[ids == id, ]
    single <- mt_dbbmm_variance(sub, location_error = sigma_by_id[[id]],
                                 window_size = 31, margin = 11)
    expect_equal(multi[[id]]$variance, single$variance,
                 tolerance = .Machine$double.eps^0.5,
                 info = paste("track", id))
  }
})

test_that("mt_dbbmm_variance stores per-fix location_error on the variance object", {
  x <- .make_fisher_multi()
  loc_err <- rep(25, nrow(x))
  res <- mt_dbbmm_variance(x, location_error = loc_err,
                            window_size = 31, margin = 11)

  for (nm in names(res)) {
    stored <- res[[nm]]$track_data$location_error
    expect_type(stored, "double")
    expect_length(stored, res[[nm]]$track_data$n_locs)
    expect_true(all(stored == 25))
  }
})

test_that("mt_dbbmm_variance rejects a wrong-length location_error vector", {
  x <- .make_fisher_multi()
  expect_error(
    mt_dbbmm_variance(x, location_error = rep(25, nrow(x) - 1L),
                       window_size = 31, margin = 11),
    "length"
  )
})

test_that("NA imputation: median default fills NAs in per-fix vector", {
  x <- .make_fisher_multi()
  ids <- as.character(mt_track_id(x))
  loc_err <- rep(25, nrow(x))
  ## Inject NAs into ~10% of rows of a single track.
  target_id <- ids[1L]
  target_rows <- which(ids == target_id)
  na_idx <- target_rows[seq(1L, length(target_rows), by = 10L)]
  loc_err[na_idx] <- NA_real_

  res <- mt_dbbmm_variance(x, location_error = loc_err,
                            window_size = 31, margin = 11)
  stored <- res[[target_id]]$track_data$location_error
  expect_false(any(is.na(stored)))
  expect_true(all(stored == 25))  # median of all-25 non-NA values is 25
})

test_that("NA imputation: explicit 'zero' fills NAs with 0", {
  x <- .make_fisher_multi()
  loc_err <- rep(25, nrow(x))
  loc_err[1:5] <- NA_real_

  res <- mt_dbbmm_variance(x, location_error = loc_err,
                            window_size = 31, margin = 11,
                            location_error_na = "zero")
  first <- res[[1L]]
  expect_equal(first$track_data$location_error[1:5], rep(0, 5))
})

test_that("UD list dispatch uses stored per-fix location_error by default", {
  x <- .make_fisher_multi()
  ids <- as.character(mt_track_id(x))
  uid <- unique(ids)
  sigma_by_id <- setNames(seq(10, 10 + 5 * (length(uid) - 1), by = 5), uid)
  loc_err <- sigma_by_id[ids]

  var_list <- mt_dbbmm_variance(x, location_error = loc_err,
                                 window_size = 31, margin = 11)
  ud <- mt_dbbmm_ud(var_list, dim_size = 50, verbose = FALSE)

  expect_s4_class(ud, "SpatRaster")
  expect_equal(terra::nlyr(ud), length(uid))
  expect_setequal(names(ud), names(var_list))
  ## Each layer should sum to ~1
  for (i in seq_len(terra::nlyr(ud))) {
    s <- sum(terra::values(ud[[i]]), na.rm = TRUE)
    expect_equal(s, 1, tolerance = 1e-6, info = paste("layer", i))
  }
})

test_that("UD list dispatch rejects per-row location_error vectors", {
  x <- .make_fisher_multi()
  var_list <- mt_dbbmm_variance(x, location_error = 25,
                                 window_size = 31, margin = 11)
  expect_error(
    mt_dbbmm_ud(var_list, location_error = rep(25, nrow(x)),
                 dim_size = 50, verbose = FALSE),
    "list-dispatch",
    fixed = FALSE
  )
})

test_that("UD list dispatch accepts a scalar location_error override", {
  x <- .make_fisher_multi()
  var_list <- mt_dbbmm_variance(x, location_error = 25,
                                 window_size = 31, margin = 11)
  ud <- suppressMessages(suppressWarnings(
    mt_dbbmm_ud(var_list, location_error = 30,
                 dim_size = 50, verbose = FALSE)
  ))
  expect_s4_class(ud, "SpatRaster")
})

## ---- dBGB twin tests ----------------------------------------------

test_that("mt_dbgb_variance accepts a per-fix vector on multi-track input", {
  x <- .make_fisher_multi()
  loc_err <- rep(25, nrow(x))
  res <- suppressMessages(suppressWarnings(
    mt_dbgb_variance(x, location_error = loc_err,
                      window_size = 31, margin = 15)
  ))
  expect_true(all(vapply(res, inherits, logical(1), "mt_dbgb_variance")))
})

test_that("mt_dbgb_variance per-track slice matches direct single-track fit", {
  x <- .make_fisher_multi()
  ids <- as.character(mt_track_id(x))
  uid <- unique(ids)
  sigma_by_id <- setNames(seq(10, 10 + 5 * (length(uid) - 1), by = 5), uid)
  loc_err <- sigma_by_id[ids]

  multi <- mt_dbgb_variance(x, location_error = loc_err,
                             window_size = 31, margin = 15)
  for (id in names(multi)) {
    sub <- x[ids == id, ]
    single <- mt_dbgb_variance(sub, location_error = sigma_by_id[[id]],
                                window_size = 31, margin = 15)
    expect_equal(multi[[id]]$para_sd, single$para_sd,
                 tolerance = .Machine$double.eps^0.5,
                 info = paste("track", id, "para_sd"))
    expect_equal(multi[[id]]$orth_sd, single$orth_sd,
                 tolerance = .Machine$double.eps^0.5,
                 info = paste("track", id, "orth_sd"))
  }
})

test_that("mt_dbgb_variance stores per-fix location_error on the variance object", {
  x <- .make_fisher_multi()
  loc_err <- rep(25, nrow(x))
  res <- mt_dbgb_variance(x, location_error = loc_err,
                           window_size = 31, margin = 15)
  for (nm in names(res)) {
    stored <- res[[nm]]$track_data$location_error
    expect_length(stored, res[[nm]]$track_data$n_locs)
    expect_true(all(stored == 25))
  }
})

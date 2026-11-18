make_fake_dbbmm <- function(n = 20) {
  structure(
    list(
      variance   = rep(1, n),
      in_windows = rep(31, n),
      interest   = rep(TRUE, n),
      break_list = integer(0),
      window_size = 31,
      margin      = 11,
      track_data  = list(x = seq_len(n), y = seq_len(n),
                         time_mins = seq_len(n), n_locs = n)
    ),
    class = "mt_dbbmm_variance"
  )
}

make_fake_dbgb <- function(n = 20) {
  structure(
    list(
      para_sd      = rep(1, n),
      orth_sd      = rep(1, n),
      n_estim      = rep(30, n),
      seg_interest = rep(TRUE, n),
      margin       = 15,
      window_size  = 31,
      track_data   = list(x = seq_len(n), y = seq_len(n),
                          time_mins = seq_len(n), n_locs = n)
    ),
    class = "mt_dbgb_variance"
  )
}

test_that("mt_mask_segments clears interest on dBBMM by integer index", {
  v <- make_fake_dbbmm(20)
  out <- mt_mask_segments(v, c(3, 7, 15))
  expect_false(out$interest[3])
  expect_false(out$interest[7])
  expect_false(out$interest[15])
  expect_true(all(out$interest[-c(3, 7, 15)]))
})

test_that("mt_mask_segments clears seg_interest on dBGB by integer index", {
  v <- make_fake_dbgb(20)
  out <- mt_mask_segments(v, c(3, 7, 15))
  expect_false(out$seg_interest[3])
  expect_false(out$seg_interest[7])
  expect_false(out$seg_interest[15])
  expect_true(all(out$seg_interest[-c(3, 7, 15)]))
})

test_that("dBBMM and dBGB masks are indexed identically", {
  vb <- make_fake_dbbmm(20)
  vd <- make_fake_dbgb(20)
  idx <- c(2L, 5L, 11L)
  expect_equal(
    mt_mask_segments(vb, idx)$interest,
    mt_mask_segments(vd, idx)$seg_interest
  )
})

test_that("logical segments of length n-1 are accepted", {
  v <- make_fake_dbbmm(10)
  mask <- c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
  out <- mt_mask_segments(v, mask)
  expect_false(out$interest[2])
  expect_false(out$interest[5])
  expect_true(all(out$interest[-c(2, 5)]))
})

test_that("logical segments of length n with trailing NA are accepted", {
  v <- make_fake_dbbmm(10)
  mask <- c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, NA)
  out <- mt_mask_segments(v, mask)
  expect_false(out$interest[2])
  expect_false(out$interest[5])
  expect_true(all(out$interest[-c(2, 5)]))
})

test_that("NA integer indices are dropped silently", {
  v <- make_fake_dbbmm(10)
  out <- mt_mask_segments(v, c(3L, NA_integer_, 7L))
  expect_false(out$interest[3])
  expect_false(out$interest[7])
  expect_equal(sum(!out$interest), 2L)
})

test_that("empty segments leave the object untouched", {
  v <- make_fake_dbbmm(10)
  expect_equal(mt_mask_segments(v, integer(0))$interest, v$interest)
  expect_equal(
    mt_mask_segments(v, rep(FALSE, 9))$interest,
    v$interest
  )
})

test_that("out-of-range indices raise an informative error", {
  v <- make_fake_dbbmm(10)
  expect_error(mt_mask_segments(v, c(1L, 15L)), "must be integers in 1:9")
  expect_error(mt_mask_segments(v, 0L),          "must be integers in 1:9")
})

test_that("logical vector with wrong length raises an informative error", {
  v <- make_fake_dbbmm(10)
  expect_error(mt_mask_segments(v, rep(FALSE, 5)), "length n-1")
})

test_that("non-numeric non-logical segments raise an error", {
  v <- make_fake_dbbmm(10)
  expect_error(mt_mask_segments(v, "foo"), "integer or logical")
})

test_that("default method errors on unknown class", {
  expect_error(
    mt_mask_segments(list(a = 1), 1L),
    "mt_dbbmm_variance.*mt_dbgb_variance"
  )
})

test_that("list method with shared segments broadcasts to each element", {
  lst <- list(A = make_fake_dbbmm(20), B = make_fake_dbbmm(25))
  out <- mt_mask_segments(lst, c(4L, 9L))
  expect_false(out$A$interest[4])
  expect_false(out$A$interest[9])
  expect_false(out$B$interest[4])
  expect_false(out$B$interest[9])
})

test_that("list method with named list of segments applies per element", {
  lst <- list(A = make_fake_dbbmm(20), B = make_fake_dbbmm(25))
  out <- mt_mask_segments(lst, list(A = 4L, B = c(10L, 11L)))
  expect_false(out$A$interest[4])
  expect_true(out$A$interest[10])
  expect_false(out$B$interest[10])
  expect_false(out$B$interest[11])
})

test_that("dBGB kernel inputs keep every fix and gate the masked bridge", {
  # Since 0.4.6 masking a dBGB bridge works as for dBBMM: every core fix
  # stays in the kernel input and only the masked bridge is skipped.
  # Before, the input was subset to seg_interest, which dropped fix i and
  # bridged i-1 to i+1 across it.  The rev-OR guard (an earlier
  # silent-undo bug) is kept: an asymmetric mask must stay asymmetric.
  ns <- asNamespace("move2utils")
  v <- make_fake_dbgb(20)
  k0 <- ns$.dbgb_kernel_inputs(v, location_error = 1, time_step = NULL)
  vm <- mt_mask_segments(v, 3L)
  k1 <- ns$.dbgb_kernel_inputs(vm, location_error = 1, time_step = NULL)
  expect_identical(k1$x, k0$x)
  expect_identical(k1$t, k0$t)
  expect_true(all(k0$interest))
  expect_identical(which(!k1$interest), 3L)
  # a core that starts later: segment indices still refer to the full track
  v$n_estim[1:2] <- NA
  v$seg_interest[1:2] <- FALSE
  vm <- mt_mask_segments(v, 5L)
  k2 <- ns$.dbgb_kernel_inputs(vm, location_error = 1, time_step = NULL)
  expect_identical(length(k2$x), 18L)
  expect_identical(which(!k2$interest), 3L)   # bridge 5 of the full track
  # every bridge masked: refused
  expect_error(ns$.dbgb_kernel_inputs(mt_mask_segments(v, 1:19),
                                      location_error = 1, time_step = NULL),
               class = "move2utils_ud_too_few_fixes")
})

test_that("masking a dBGB bridge removes exactly that bridge from the UD", {
  skip_on_cran()
  skip_if_not_installed("move2")
  ns <- asNamespace("move2utils")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  leroy <- move2::filter_track_data(fishers, .track_id = "M1")
  leroy <- sf::st_transform(leroy, move2::mt_aeqd_crs(leroy))
  v <- mt_dbgb_variance(leroy, window_size = 31, margin = 11,
                        location_error = 20)
  core <- which(v$seg_interest)
  lag <- diff(v$track_data$time_mins)
  ## the longest gap inside the core
  i <- core[-length(core)][which.max(lag[core[-length(core)]])]
  j <- i - core[1] + 1L                        # its index in the core
  g <- mt_dbgb_ud(v, location_error = 20, res = 200, verbose = FALSE)
  g <- terra::extend(g, 20)
  gv <- ns$.ud_grid_vectors(g, "dbgb")
  kin <- ns$.dbgb_kernel_inputs(v, 20, NULL)
  raw <- function(intr) {
    .Call("bgb_omp", kin$x, kin$y, kin$s1, kin$s2, kin$t, kin$le,
          gv$x, gv$y, kin$time_step, kin$ext, intr, PACKAGE = "move2utils")
  }
  only <- rep(FALSE, length(kin$x)); only[j] <- TRUE
  all_b <- raw(kin$interest)
  ## through the public path: same fixes, bridge j off
  kin_m <- ns$.dbgb_kernel_inputs(mt_mask_segments(v, i), 20, NULL)
  expect_identical(kin_m$x, kin$x)
  minus <- raw(kin_m$interest)
  expect_equal(minus + raw(only), all_b, tolerance = 1e-10)
  ## and the UD loses mass where the gap's bridge ran
  u0 <- mt_dbgb_ud(v, raster = g, location_error = 20, verbose = FALSE)
  u1 <- mt_dbgb_ud(mt_mask_segments(v, i), raster = g, location_error = 20,
                   verbose = FALSE)
  td <- v$track_data
  mid <- cbind((td$x[i] + td$x[i + 1]) / 2, (td$y[i] + td$y[i + 1]) / 2)
  expect_lt(terra::extract(u1, mid)[1, 1], terra::extract(u0, mid)[1, 1])
})

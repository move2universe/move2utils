## Cross-package numerical parity: dBBMM variance.
##
## Compares mt_dbbmm_variance() output against the legacy
## move::brownian.motion.variance.dyn() reference on the same input
## (move::leroy -- the canonical fisher track that both packages can
## consume).  This is a permanent regression guard against drift between
## the move2utils C kernel (src/bm_variance_c.c) and the move-era kernel
## it was migrated from.
##
## Tolerance: the per-fix variance estimate agrees with move at
## machine precision (median relative ~1e-8) for the vast majority of
## fixes.  A subset of stationary-segment fixes hit move's lower-bound
## clamp at ~5.5e-5 while our kernel converges to ~1e-10; both
## implementations agree the variance is effectively zero, they differ
## on the implementation floor.  We assert on the median relative diff
## rather than the max, which correctly captures the "C kernel == move
## kernel at machine precision" story without being fooled by the
## handful of clamp-floor mismatches.  See
## benchmarks/c_port_vs_move/REPORT.md for the broader numerical-
## agreement study on the synthetic CPF tracks.

test_that("mt_dbbmm_variance matches move::brownian.motion.variance.dyn on move::leroy", {
  skip_on_cran()
  skip_if_not_installed("move")
  skip_if_not_installed("move2")
  skip_if_not_installed("sp")

  ## Fetch leroy from the move package and coerce to move2.
  e <- new.env()
  utils::data("leroy", package = "move", envir = e)
  leroy_move <- e$leroy
  leroy_m2   <- move2::mt_as_move2(leroy_move)

  ## Project both to the same AEQD so the two pipelines see the same
  ## projected coordinates and any disagreement is attributable to the
  ## estimator alone.
  aeqd_crs <- move2::mt_aeqd_crs(leroy_m2)
  leroy_m2_proj <- sf::st_transform(leroy_m2, aeqd_crs)
  cc <- sf::st_coordinates(leroy_m2_proj)

  ## Re-build the move::Move from the projected coordinates so the
  ## legacy pipeline sees the *same* numerical input.
  leroy_move_proj <- move::move(
    x      = cc[, 1],
    y      = cc[, 2],
    time   = move2::mt_time(leroy_m2_proj),
    proj   = sf::st_crs(leroy_m2_proj)$proj4string,
    animal = "leroy"
  )

  ws     <- 31L
  margin <- 11L
  loc_err <- 25

  ref_obj <- move::brownian.motion.variance.dyn(
    leroy_move_proj,
    location.error = loc_err, window.size = ws, margin = margin)
  ref_var <- ref_obj@means

  our_obj <- mt_dbbmm_variance(
    leroy_m2_proj,
    location_error = loc_err, window_size = ws, margin = margin)
  our_var <- mt_motion_variance(our_obj)
  if (is.data.frame(our_var)) our_var <- our_var[[1]]

  ok <- !is.na(ref_var) & !is.na(our_var)
  expect_equal(is.na(ref_var), is.na(our_var))

  abs_diff <- abs(ref_var[ok] - our_var[ok])
  rel_diff <- abs_diff / pmax(abs(ref_var[ok]), 1e-12)
  message(sprintf(
    "dBBMM variance vs move::leroy: n_ok=%d  med|d|=%.3e  max|d|=%.3e  med rel=%.3e  q90 rel=%.3e",
    sum(ok), median(abs_diff), max(abs_diff), median(rel_diff), quantile(rel_diff, .90)))

  ## Median relative agreement at machine precision: the same Brent's-
  ## method kernel converges to the same minimum on the vast majority
  ## of fixes.
  expect_lt(median(rel_diff), 1e-6)
  ## Maximum absolute variance difference stays small even where the
  ## two implementations clamp differently on zero-velocity segments.
  expect_lt(max(abs_diff), 1e-3)
})

## Cross-package numerical parity: dBGB variance.
##
## Compares mt_dbgb_variance() output against the legacy
## move::dynBGBvariance() reference on move::leroy.
##
## Tolerance scope: optimiser-convergence floor, not machine precision.
## The dBGB single-window estimator does a per-axis Brent's-method
## optimisation against a non-convex surface with an inner breakpoint
## search; the legacy implementation does the equivalent in pure R via
## optim() + a different breakpoint discipline.  Both converge to the
## same minimum modulo their internal tolerances (~1e-3 on sigma for
## well-behaved windows; up to a few % relative at fixes inside
## degenerate-variance regions where sigma is enormous and the
## optimisation surface is flat).  We assert on the median absolute
## difference -- a single-fix outlier driven by optimiser-floor noise
## should not break the regression.  See
## benchmarks/c_port_vs_move/REPORT.md for the broader analysis.

test_that("mt_dbgb_variance matches move::dynBGBvariance on move::leroy (median floor)", {
  skip_on_cran()
  skip_if_not_installed("move")
  skip_if_not_installed("move2")
  skip_if_not_installed("sp")

  e <- new.env()
  utils::data("leroy", package = "move", envir = e)
  leroy_move <- e$leroy
  leroy_m2   <- move2::mt_as_move2(leroy_move)

  aeqd_crs <- move2::mt_aeqd_crs(leroy_m2)
  leroy_m2_proj <- sf::st_transform(leroy_m2, aeqd_crs)
  cc <- sf::st_coordinates(leroy_m2_proj)

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

  ref_obj <- suppressWarnings(move::dynBGBvariance(
    leroy_move_proj, locErr = loc_err, margin = margin, windowSize = ws))
  ref_para <- ref_obj@paraSd
  ref_orth <- ref_obj@orthSd

  our_obj <- mt_dbgb_variance(
    leroy_m2_proj,
    location_error = loc_err, margin = margin, window_size = ws)
  mv <- mt_motion_variance(our_obj)
  ## move's slots are SIGMA; our object stores VARIANCE -> take sqrt.
  our_para <- sqrt(mv$para)
  our_orth <- sqrt(mv$orth)

  ok_p <- !is.na(ref_para) & !is.na(our_para)
  ok_o <- !is.na(ref_orth) & !is.na(our_orth)
  expect_equal(is.na(ref_para), is.na(our_para))
  expect_equal(is.na(ref_orth), is.na(our_orth))

  d_p <- abs(ref_para[ok_p] - our_para[ok_p])
  d_o <- abs(ref_orth[ok_o] - our_orth[ok_o])
  message(sprintf(
    "dBGB sigma_para vs move::leroy: n_ok=%d  med|d|=%.3e  q90|d|=%.3e  q99|d|=%.3e  max=%.3e",
    sum(ok_p), median(d_p), quantile(d_p, .90), quantile(d_p, .99), max(d_p)))
  message(sprintf(
    "dBGB sigma_orth vs move::leroy: n_ok=%d  med|d|=%.3e  q90|d|=%.3e  q99|d|=%.3e  max=%.3e",
    sum(ok_o), median(d_o), quantile(d_o, .90), quantile(d_o, .99), max(d_o)))

  ## Bulk agreement: 90 % of fixes within 1e-2 on sigma (~ centimetre).
  expect_lt(quantile(d_p, .90), 1e-2)
  expect_lt(quantile(d_o, .90), 1e-2)
  ## Median agreement at optimiser-tolerance level.
  expect_lt(median(d_p), 1e-3)
  expect_lt(median(d_o), 1e-3)
})

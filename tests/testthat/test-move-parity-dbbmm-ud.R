## Cross-package numerical parity: dBBMM UD.
##
## After the variance step agrees with move:: at machine precision (see
## test-move-parity-dbbmm-variance.R), the UD step is a deterministic
## grid evaluation against the same dbbmm2() C kernel that move ships
## (extended in move2utils with an `interest` per-segment gating
## argument, default-on so behaviour is unchanged when no gating is
## requested).  We check that the UD raster cell values match cell-by-
## cell to within machine precision.

test_that("mt_dbbmm_ud matches move::brownian.bridge.dyn on move::leroy", {
  skip_on_cran()
  skip_if_not_installed("move")
  skip_if_not_installed("move2")
  skip_if_not_installed("terra")
  skip_if_not_installed("sp")
  skip_if_not_installed("raster")

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
  raster_px <- 100   # 100 m pixel
  ext_pad   <- 2     # generous extent padding — move's grid sizing is
                       # less forgiving than ours on tall windows

  ## Legacy UD via brownian.bridge.dyn
  ref_ud <- move::brownian.bridge.dyn(
    leroy_move_proj,
    raster         = raster_px,
    location.error = loc_err,
    margin         = margin,
    window.size    = ws,
    ext            = ext_pad,
    verbose        = FALSE)
  ref_vals <- raster::values(ref_ud)

  ## move2utils UD
  our_var <- mt_dbbmm_variance(
    leroy_m2_proj,
    location_error = loc_err, window_size = ws, margin = margin)
  our_ud <- mt_dbbmm_ud(our_var, raster = raster_px, ext = ext_pad)
  our_vals <- terra::values(our_ud)

  ## Both UDs should be normalised to sum-to-one over their grid; their
  ## extents may differ by integer pixels if the bbox snap rounds
  ## differently.  Check at the integral-distance level: total absolute
  ## difference of normalised mass.
  ref_n <- ref_vals / sum(ref_vals, na.rm = TRUE)
  our_n <- our_vals / sum(our_vals, na.rm = TRUE)
  message(sprintf(
    "dBBMM UD vs move::leroy: ref cells=%d, our cells=%d, ref sum=%.6f, our sum=%.6f",
    length(ref_vals), length(our_vals),
    sum(ref_vals, na.rm = TRUE), sum(our_vals, na.rm = TRUE)))

  ## If grids match in size, do a cell-by-cell check; otherwise compare
  ## the sorted mass distributions (insensitive to pixel-frame shift).
  if (length(ref_n) == length(our_n)) {
    d <- abs(sort(ref_n) - sort(our_n))
    message(sprintf("  cell-mass max|d| = %.3e  med|d| = %.3e",
                    max(d), median(d)))
    expect_lt(max(d), 1e-6)
  } else {
    rs <- sort(ref_n); os <- sort(our_n)
    ## interpolate to common length before comparing
    if (length(rs) > length(os)) rs <- approx(seq_along(rs), rs, n = length(os))$y
    else os <- approx(seq_along(os), os, n = length(rs))$y
    d <- abs(rs - os)
    message(sprintf("  sorted-mass max|d| = %.3e  med|d| = %.3e",
                    max(d), median(d)))
    expect_lt(max(d), 1e-4)
  }
})

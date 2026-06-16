## Cross-package numerical parity: Earth-mover's distance.
##
## move2utils::emd(method = "exact") delegates to emdist::emd() with
## the same simplex-LP formulation that move::emd() uses under the
## hood.  On the canonical move-package example (the dbbmmstack UD
## stack) the two implementations agree on the scalar transport
## distance at the LP solver's convergence floor (~1e-4 relative); the
## residual is whichever solver instance hit its internal iteration
## cap first.  We assert agreement to 1e-3 (0.1 %), comfortably above
## the measured floor.

test_that("emd(method='exact') matches move::emd on dbbmmstack", {
  skip_on_cran()
  skip_if_not_installed("move")
  skip_if_not_installed("emdist")
  skip_if_not_installed("terra")
  skip_if_not_installed("raster")

  e <- new.env()
  utils::data("dbbmmstack", package = "move", envir = e)
  stk <- e$dbbmmstack
  ## dbbmmstack ships as a RasterStack of two UDs.
  expect_true(raster::nlayers(stk) >= 2L)

  ## --- legacy move::emd on the original RasterStack ---
  d_move <- as.numeric(move::emd(stk))

  ## --- move2utils::emd(method="exact") on the same data as terra ---
  ## Rebuild the stack as a SpatRaster (terra::rast() can't coerce
  ## move's DBBMMStack class directly through raster's CRS attribute).
  ext_v <- as.vector(raster::extent(stk))
  stk_terra <- terra::rast(
    nrows = nrow(stk), ncols = ncol(stk),
    xmin = ext_v[1], xmax = ext_v[2], ymin = ext_v[3], ymax = ext_v[4],
    crs  = as.character(raster::crs(stk)),
    nlyr = raster::nlayers(stk))
  for (i in seq_len(raster::nlayers(stk))) {
    terra::values(stk_terra[[i]]) <- raster::values(stk)[, i]
  }
  names(stk_terra) <- names(stk)
  d_ours <- as.numeric(emd(stk_terra, method = "exact", mask_quantile = 1))

  message(sprintf("emd exact: move = %.6e   move2utils = %.6e   diff = %.3e",
                  d_move, d_ours, abs(d_move - d_ours)))

  ## LP-vs-LP at the solver's convergence floor.
  expect_lt(abs(d_move - d_ours) / max(abs(d_move), 1e-12), 1e-3)
})

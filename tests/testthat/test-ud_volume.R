test_that("ud_volume rejects non-SpatRaster input", {
  expect_error(ud_volume(matrix(1:4, 2, 2)),
               class = "move2utils_ud_volume_not_spatraster")
})

test_that("ud_volume outputs values in [0, 1] with max ~1", {
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 1,
                   ymin = 0, ymax = 1, vals = runif(400))
  r <- r / sum(terra::values(r))        # normalise to a UD
  v <- ud_volume(r)

  vals <- terra::values(v)[, 1]
  expect_true(all(vals >= 0 & vals <= 1 + 1e-9, na.rm = TRUE))
  expect_equal(max(vals, na.rm = TRUE), 1, tolerance = 1e-9)
})

test_that("largest UD cell receives the smallest volume quantile", {
  set.seed(1)
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1,
                   ymin = 0, ymax = 1, vals = runif(100))
  r <- r / sum(terra::values(r))

  v <- ud_volume(r)
  i_max <- which.max(terra::values(r)[, 1])
  vols  <- terra::values(v)[, 1]

  expect_equal(vols[i_max], min(vols, na.rm = TRUE))
})

test_that("ud_volume preserves NA cells", {
  r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 1,
                   ymin = 0, ymax = 1, vals = runif(25))
  r <- r / sum(terra::values(r))
  terra::values(r)[c(1, 5, 13)] <- NA

  v <- ud_volume(r)
  expect_equal(which(is.na(terra::values(v)[, 1])),
               c(1, 5, 13))
})

test_that("ud_volume handles multi-layer input layer-by-layer", {
  set.seed(2)
  r1 <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1,
                    ymin = 0, ymax = 1, vals = runif(100))
  r1 <- r1 / sum(terra::values(r1))
  r2 <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1,
                    ymin = 0, ymax = 1, vals = runif(100))
  r2 <- r2 / sum(terra::values(r2))
  stk <- c(r1, r2)
  names(stk) <- c("a", "b")

  v <- ud_volume(stk)
  expect_equal(terra::nlyr(v), 2L)
  expect_equal(names(v), c("a", "b"))
  for (i in seq_len(terra::nlyr(v))) {
    expect_equal(max(terra::values(v[[i]])[, 1], na.rm = TRUE), 1,
                 tolerance = 1e-9)
  }
})

test_that("ud_volume matches the legacy getVolumeUD inline recipe", {
  set.seed(3)
  r <- terra::rast(nrows = 15, ncols = 15, xmin = 0, xmax = 1,
                   ymin = 0, ymax = 1, vals = runif(225))
  r <- r / sum(terra::values(r))

  ## inline recipe used in the vignettes before ud_volume()
  vol_inline <- function(ud) {
    v <- terra::values(ud)[, 1]
    o <- order(v, decreasing = TRUE)
    v[o] <- cumsum(v[o])
    out <- ud
    terra::values(out) <- v
    out
  }
  v_fun    <- ud_volume(r)
  v_inline <- vol_inline(r)

  expect_equal(terra::values(v_fun)[, 1],
               terra::values(v_inline)[, 1],
               tolerance = 1e-12)
})

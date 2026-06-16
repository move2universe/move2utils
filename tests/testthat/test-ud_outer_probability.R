test_that("ud_outer_probability validates inputs", {
  r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 1,
                   ymin = 0, ymax = 1, vals = runif(25))
  r <- r / sum(terra::values(r))
  expect_error(ud_outer_probability("nope", r),
               "sf / sfc / move2")
  expect_error(ud_outer_probability(matrix(1:5, ncol = 1), r),
               "at least two columns")
  expect_error(ud_outer_probability(data.frame(x = 0.5, y = 0.5),
                                     "not a raster"),
               "SpatRaster")
})

test_that("ud_outer_probability returns values in [0, 1]", {
  set.seed(1)
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 10,
                   ymin = 0, ymax = 10, vals = runif(400))
  r <- r / sum(terra::values(r))
  pts <- cbind(runif(50, 0, 10), runif(50, 0, 10))
  op <- ud_outer_probability(pts, r)
  expect_true(all(op >= 0 & op <= 1, na.rm = TRUE))
})

test_that("peak-density location has the smallest outer probability", {
  ## hand-crafted UD with a clear peak at the centre
  grid <- expand.grid(x = seq(0.05, 0.95, by = 0.1),
                      y = seq(0.05, 0.95, by = 0.1))
  vals <- exp(-((grid$x - 0.5)^2 + (grid$y - 0.5)^2) * 20)
  vals <- vals / sum(vals)
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1,
                   ymin = 0, ymax = 1, vals = vals)

  pts <- rbind(
    c(0.5, 0.5),   # peak
    c(0.05, 0.05)  # corner
  )
  op <- ud_outer_probability(pts, r)
  ## peak has smaller outer probability than corner (core area)
  expect_lt(op[1], op[2])
  expect_lt(op[1], 0.5)
  expect_gt(op[2], 0.5)
})

test_that("ud_outer_probability accepts a move2 object", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  leroy <- fishers[move2::mt_track_id(fishers) == "M4", ][1:100, ]
  leroy <- sf::st_transform(leroy, move2::mt_aeqd_crs(leroy))

  ud <- terra::rast(nrows = 10, ncols = 10,
                    xmin = min(sf::st_coordinates(leroy)[, 1]) - 1000,
                    xmax = max(sf::st_coordinates(leroy)[, 1]) + 1000,
                    ymin = min(sf::st_coordinates(leroy)[, 2]) - 1000,
                    ymax = max(sf::st_coordinates(leroy)[, 2]) + 1000,
                    vals = runif(100),
                    crs = sf::st_crs(leroy)$wkt)
  ud <- ud / sum(terra::values(ud))

  op <- ud_outer_probability(leroy, ud)
  expect_length(op, nrow(leroy))
  expect_true(all(op >= 0 & op <= 1, na.rm = TRUE))
})

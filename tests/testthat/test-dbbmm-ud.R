test_that("mt_dbbmm_ud returns a SpatRaster", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  leroy <- fishers[mt_track_id(fishers) == "F1", ]
  bb <- st_bbox(leroy)
  leroy_proj <- st_transform(leroy, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  # Use small window for speed, sufficient extent for grid
  ud <- mt_dbbmm_ud(leroy_proj, location_error = 25,
                  window_size = 31, margin = 11,
                  ext = 1.0, dim_size = 100, verbose = FALSE)

  expect_s4_class(ud, "SpatRaster")
  expect_true(abs(sum(terra::values(ud), na.rm = TRUE) - 1) < 0.01)
})

test_that("mt_dbbmm_ud accepts pre-computed variance object", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  leroy <- fishers[mt_track_id(fishers) == "F1", ]
  bb <- st_bbox(leroy)
  leroy_proj <- st_transform(leroy, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  var_obj <- mt_dbbmm_variance(leroy_proj, location_error = 25,
                                 window_size = 31, margin = 11)
  ud <- mt_dbbmm_ud(var_obj, location_error = 25, ext = 1.0,
                  dim_size = 100, verbose = FALSE)

  expect_s4_class(ud, "SpatRaster")
  expect_true(abs(sum(terra::values(ud), na.rm = TRUE) - 1) < 0.01)
})

test_that("mt_dbbmm_ud accepts numeric cell size", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  leroy <- fishers[mt_track_id(fishers) == "F1", ]
  bb <- st_bbox(leroy)
  leroy_proj <- st_transform(leroy, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  ud <- mt_dbbmm_ud(leroy_proj, raster = 100, location_error = 25,
                  window_size = 31, margin = 11, ext = 1.0, verbose = FALSE)

  expect_s4_class(ud, "SpatRaster")
  # Cell size should be approximately 100m
  expect_equal(terra::res(ud)[1], 100, tolerance = 1)
})

test_that("multi-track input returns a named list", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  bb <- st_bbox(fishers)
  fishers_proj <- st_transform(fishers, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  # Should return a list, not error
  result <- mt_dbbmm_variance(fishers_proj, location_error = 25,
                               window_size = 31, margin = 11)
  expect_type(result, "list")
  # Each element should be an mt_dbbmm_variance
  expect_true(all(vapply(result, inherits, logical(1), "mt_dbbmm_variance")))
  # Named by track ID
  expect_true(all(names(result) %in% unique(as.character(mt_track_id(fishers_proj)))))
})

test_that("multi-track UD returns multi-layer SpatRaster", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  # Use just 2 tracks for speed
  two <- fishers[mt_track_id(fishers) %in% c("F1", "F2"), ]
  bb <- st_bbox(two)
  two_proj <- st_transform(two, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  var_list <- mt_dbbmm_variance(two_proj, location_error = 25,
                                 window_size = 31, margin = 11)
  ud_stack <- mt_dbbmm_ud(var_list, location_error = 25,
                            ext = 1.0, dim_size = 100, verbose = FALSE)

  expect_s4_class(ud_stack, "SpatRaster")
  expect_equal(terra::nlyr(ud_stack), 2)
  expect_true(all(names(ud_stack) %in% c("F1", "F2")))
  # Each layer should sum to 1
  for (i in 1:terra::nlyr(ud_stack)) {
    expect_true(abs(sum(terra::values(ud_stack[[i]]), na.rm = TRUE) - 1) < 0.01)
  }
})

test_that("tracks with too few locations are skipped with warning", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  # Take 2 tracks, one with plenty of data and one with very few
  f1 <- fishers[mt_track_id(fishers) == "F1", ]
  f2 <- fishers[mt_track_id(fishers) == "F2", ][1:10, ]  # only 10 locations
  combined <- rbind(f1, f2)
  bb <- st_bbox(combined)
  combined_proj <- st_transform(combined, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  expect_warning(
    result <- mt_dbbmm_variance(combined_proj, location_error = 25,
                                 window_size = 31, margin = 11),
    "skipped"
  )
  # Only F1 should be in the results
  expect_length(result, 1)
  expect_true("F1" %in% names(result))
})

test_that("empty geometries are rejected", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  f1 <- fishers[mt_track_id(fishers) == "F1", ]
  bb <- st_bbox(f1[!st_is_empty(f1), ])
  f1_proj <- st_transform(f1, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  expect_error(
    mt_dbbmm_variance(f1_proj, location_error = 25,
                       window_size = 31, margin = 11),
    "empty geometries"
  )
})

test_that("lon/lat input is auto-projected", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  f1 <- fishers[mt_track_id(fishers) == "F1", ]

  ## lon/lat is reprojected to a local AEQD internally rather than rejected
  ## (projection stratification 2026-06-08; was an error matching
  ## "projected CRS").
  expect_message(
    v <- mt_dbbmm_variance(f1, location_error = 25,
                       window_size = 31, margin = 11),
    "Auto-projecting"
  )
  expect_s3_class(v, "mt_dbbmm_variance")
})

test_that("window_size < 2*margin+1 is rejected", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  f1 <- fishers[mt_track_id(fishers) == "F1", ]
  bb <- st_bbox(f1)
  f1_proj <- st_transform(f1, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  expect_error(
    mt_dbbmm_variance(f1_proj, location_error = 25,
                       window_size = 9, margin = 5),
    "2 \\* margin"
  )
})

test_that("non-move2 input is rejected", {
  expect_error(
    mt_dbbmm_variance(data.frame(x = 1:10), location_error = 25,
                       window_size = 5, margin = 3),
    "move2 object"
  )
})

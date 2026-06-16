test_that("mt_dbgb_variance returns expected structure", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  leroy <- fishers[mt_track_id(fishers) == "F1", ]
  bb <- st_bbox(leroy)
  leroy_proj <- st_transform(leroy, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  result <- mt_dbgb_variance(leroy_proj, location_error = 25,
                               margin = 15, window_size = 31)

  expect_s3_class(result, "mt_dbgb_variance")
  expect_length(result$para_sd, nrow(leroy_proj))
  expect_length(result$orth_sd, nrow(leroy_proj))
  expect_true(any(!is.na(result$para_sd)))
  expect_true(any(!is.na(result$orth_sd)))
  expect_true(all(result$para_sd[!is.na(result$para_sd)] >= 0))
  expect_true(all(result$orth_sd[!is.na(result$orth_sd)] >= 0))
})

test_that("mt_motion_variance works on dbgb_var", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  leroy <- fishers[mt_track_id(fishers) == "F1", ]
  bb <- st_bbox(leroy)
  leroy_proj <- st_transform(leroy, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  var_obj <- mt_dbgb_variance(leroy_proj, location_error = 25,
                                margin = 15, window_size = 31)
  mv <- mt_motion_variance(var_obj)
  expect_true(is.data.frame(mv))
  expect_named(mv, c("para", "orth"))
  expect_equal(mv$para, var_obj$para_sd^2)
  expect_equal(mv$orth, var_obj$orth_sd^2)
})

test_that("mt_dbgb_ud returns a SpatRaster that sums to 1", {
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

  ud <- mt_dbgb_ud(leroy_proj, location_error = 25,
                 margin = 15, window_size = 31,
                 ext = 2.0, dim_size = 100)

  expect_s4_class(ud, "SpatRaster")
  expect_true(abs(sum(terra::values(ud), na.rm = TRUE) - 1) < 0.01)
})

test_that("mt_dbgb_variance matches move package output", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")

  fixture_file <- test_path("fixtures", "ref_dbgb_var.rds")
  skip_if_not(file.exists(fixture_file),
              "Reference fixtures not generated. Run helper-generate-fixtures.R first.")

  library(move2)
  library(sf)

  ref <- readRDS(fixture_file)
  leroy_proj <- readRDS(test_path("fixtures", "leroy_projected.rds"))

  result <- mt_dbgb_variance(leroy_proj,
                               location_error = ref$location_error,
                               margin = ref$margin,
                               window_size = ref$window_size)

  mv <- mt_motion_variance(result)
  valid_para <- !is.na(ref$para_var) & !is.na(mv$para)
  valid_orth <- !is.na(ref$orth_var) & !is.na(mv$orth)

  # Tolerance accounts for projection centering differences
  expect_equal(mv$para[valid_para], ref$para_var[valid_para],
               tolerance = 0.01)
  expect_equal(mv$orth[valid_orth], ref$orth_var[valid_orth],
               tolerance = 0.01)
})

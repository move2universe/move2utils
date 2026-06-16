test_that("mt_dbbmm_variance auto-projects lon/lat input", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  fishers <- fishers[mt_track_id(fishers) == "F1", ]

  ## A move2 object carries its CRS, so lon/lat input is reprojected to a
  ## local AEQD internally rather than rejected (projection stratification
  ## 2026-06-08; was an error matching "projected CRS").  An informative
  ## message is emitted and metric variances are produced.
  expect_message(
    v <- mt_dbbmm_variance(fishers, location_error = 25,
                             window_size = 31, margin = 11),
    "Auto-projecting"
  )
  expect_s3_class(v, "mt_dbbmm_variance")
  expect_true(length(mt_motion_variance(v)) > 0)
})

test_that("mt_dbbmm_variance works on projected data", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  leroy <- fishers[mt_track_id(fishers) == "F1", ]

  # Project to AEQD centred on track
  bb <- st_bbox(leroy)
  crs_aeqd <- st_crs(paste0("+proj=aeqd +lon_0=",
                              (bb["xmin"] + bb["xmax"]) / 2,
                              " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2,
                              " +units=m"))
  leroy_proj <- st_transform(leroy, crs_aeqd)

  result <- mt_dbbmm_variance(leroy_proj, location_error = 25,
                                window_size = 31, margin = 11)

  expect_s3_class(result, "mt_dbbmm_variance")
  expect_length(result$variance, nrow(leroy_proj))
  expect_true(any(!is.na(result$variance)))
  expect_true(all(result$variance[!is.na(result$variance)] >= 0))
  expect_true(any(result$interest))
  expect_equal(result$window_size, 31)
  expect_equal(result$margin, 11)
})

test_that("mt_dbbmm_variance errors on even window_size", {
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

  expect_error(
    mt_dbbmm_variance(leroy_proj, location_error = 25,
                        window_size = 30, margin = 11),
    "must both be odd"
  )
})

test_that("mt_dbbmm_variance errors when window larger than track", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  # Take only a few points
  tiny <- fishers[1:10, ]
  bb <- st_bbox(tiny)
  tiny_proj <- st_transform(tiny, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  expect_error(
    mt_dbbmm_variance(tiny_proj, location_error = 25,
                        window_size = 31, margin = 11),
    "cannot be larger"
  )
})

test_that("mt_motion_variance works on dbbmm_var", {
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

  var_obj <- mt_dbbmm_variance(leroy_proj, location_error = 25,
                                 window_size = 31, margin = 11)
  mv <- mt_motion_variance(var_obj)
  expect_identical(mv, var_obj$variance)
})

test_that("mt_dbbmm_variance matches move package output", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")

  fixture_file <- test_path("fixtures", "ref_dbbmm_var.rds")
  skip_if_not(file.exists(fixture_file),
              "Reference fixtures not generated. Run helper-generate-fixtures.R first.")

  library(move2)
  library(sf)

  ref <- readRDS(fixture_file)
  leroy_proj <- readRDS(test_path("fixtures", "leroy_projected.rds"))

  result <- mt_dbbmm_variance(leroy_proj,
                                location_error = ref$location_error,
                                window_size = ref$window_size,
                                margin = ref$margin)

  # Compare variance vectors — close but not identical due to
  # small differences in projection centering between sp::spTransform(center=T)
  # and sf::st_transform with manual AEQD
  valid <- !is.na(ref$variance) & !is.na(result$variance)
  # Use relative tolerance: values should agree within ~1%
  expect_equal(result$variance[valid], ref$variance[valid],
               tolerance = 0.01)
})

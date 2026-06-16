# CRS stratification for the legacy UD / variance layer (2026-06-08).
# Contract: the UD is returned in the CRS the caller supplied (long/lat in ->
# long/lat out, UTM in -> UTM out, AEQD in -> AEQD out).  The bridge math runs
# in a metric "compute" CRS: a projected metric target is used directly
# (lossless); long/lat is computed in a local AEQD and warped back + renormalised.

skip_if_not_installed("move2")
skip_if_not_installed("sf")
skip_if_not_installed("terra")

library(move2)
library(sf)

read_f1 <- function(n = 200L) {
  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  fishers[mt_track_id(fishers) == "F1", ][seq_len(n), ]
}

test_that(".crs_is_suitable accepts projected metres, rejects long/lat", {
  f1 <- read_f1()
  expect_false(.crs_is_suitable(sf::st_crs(f1)))            # long/lat
  expect_true(.crs_is_suitable(sf::st_crs(32618)))          # UTM 18N (metres)
  expect_true(.crs_is_suitable(move2::mt_aeqd_crs(f1)))     # local AEQD (metres)
})

test_that("dbbmm UD returns the caller's CRS and sums to 1 (long/lat, UTM, AEQD)", {
  f1  <- read_f1()
  utm <- sf::st_crs(32618)
  aq  <- move2::mt_aeqd_crs(f1)

  ud_ll <- suppressMessages(mt_dbbmm_ud(f1, location_error = 25,
                                        dim_size = 80, ext = 1.5))
  ud_u  <- suppressMessages(mt_dbbmm_ud(sf::st_transform(f1, utm),
                                        location_error = 25, dim_size = 80, ext = 1.5))
  ud_a  <- suppressMessages(mt_dbbmm_ud(sf::st_transform(f1, aq),
                                        location_error = 25, dim_size = 80, ext = 1.5))

  expect_true(sf::st_is_longlat(terra::crs(ud_ll)))
  expect_true(sf::st_crs(terra::crs(ud_u)) == utm)
  expect_true(sf::st_crs(terra::crs(ud_a)) == sf::st_crs(aq))
  for (u in list(ud_ll, ud_u, ud_a)) {
    expect_equal(sum(terra::values(u), na.rm = TRUE), 1, tolerance = 1e-6)
  }
})

test_that("dbbmm variance values are invariant to input CRS (long/lat vs AEQD)", {
  f1 <- read_f1()
  v_ll <- suppressMessages(mt_motion_variance(
    mt_dbbmm_variance(f1, location_error = 25, window_size = 31, margin = 11)))
  v_a  <- suppressMessages(mt_motion_variance(
    mt_dbbmm_variance(sf::st_transform(f1, move2::mt_aeqd_crs(f1)),
                      location_error = 25, window_size = 31, margin = 11)))
  ## long/lat auto-projects to the same local AEQD, so the variances match
  expect_equal(v_ll, v_a, tolerance = 1e-3)
})

test_that("multi-track dbbmm UD warps the whole stack to long/lat, each layer sums to 1", {
  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  two <- fishers[mt_track_id(fishers) %in% c("F1", "M1"), ]

  stk <- suppressMessages(mt_dbbmm_ud(two, location_error = 25,
                                      dim_size = 70, ext = 1.5))
  expect_equal(terra::nlyr(stk), 2L)
  expect_true(sf::st_is_longlat(terra::crs(stk)))
  sums <- terra::global(stk, "sum", na.rm = TRUE)[, 1]
  expect_equal(unname(sums), c(1, 1), tolerance = 1e-6)
})

test_that("an environmental template raster receives the UD on its exact grid", {
  f1  <- read_f1()
  utm <- sf::st_crs(32618)
  f1u <- sf::st_transform(f1, utm)
  bb  <- sf::st_bbox(f1u)
  pad <- 3000
  templ <- terra::rast(
    terra::ext(bb["xmin"] - pad, bb["xmax"] + pad,
               bb["ymin"] - pad, bb["ymax"] + pad),
    resolution = 120, crs = utm$wkt)

  ## metric track + metric template -> compute on the template directly
  ud_t <- suppressMessages(mt_dbbmm_ud(f1u, location_error = 25, raster = templ))
  expect_true(terra::compareGeom(ud_t, templ, stopOnError = FALSE))
  expect_equal(sum(terra::values(ud_t), na.rm = TRUE), 1, tolerance = 1e-6)

  ## long/lat track + metric template -> compute in the template CRS, return on it
  ud_x <- suppressMessages(mt_dbbmm_ud(f1, location_error = 25, raster = templ))
  expect_true(terra::compareGeom(ud_x, templ, stopOnError = FALSE))
  expect_equal(sum(terra::values(ud_x), na.rm = TRUE), 1, tolerance = 1e-6)
})

test_that("dbgb UD returns the caller's CRS and sums to 1", {
  f1 <- read_f1()
  ud_ll <- suppressMessages(mt_dbgb_ud(f1, location_error = 25,
                                       dim_size = 70, ext = 1.5))
  expect_true(sf::st_is_longlat(terra::crs(ud_ll)))
  expect_equal(sum(terra::values(ud_ll), na.rm = TRUE), 1, tolerance = 1e-6)

  utm <- sf::st_crs(32618)
  ud_u <- suppressMessages(mt_dbgb_ud(sf::st_transform(f1, utm),
                                      location_error = 25, dim_size = 70, ext = 1.5))
  expect_true(sf::st_crs(terra::crs(ud_u)) == utm)
  expect_equal(sum(terra::values(ud_u), na.rm = TRUE), 1, tolerance = 1e-6)
})

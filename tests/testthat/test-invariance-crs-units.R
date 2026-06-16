# Metamorphic / invariance guards (2026-06-08).
#
# These assert the property that a fix's scientific result must NOT depend on
# the coordinate reference system or linear unit the caller happens to supply.
# This is the bug class that slipped past every per-function correctness audit
# (the projection leak fixed in v0.4.1, and ud_outer_probability / mt_thin_distance
# found by the conceptual-integrity audit): a function can be perfectly correct
# for its input yet give a *different* correct-looking answer in another CRS.
#
# The discipline: feed the SAME track in two equivalent representations and
# assert equal output. See CODE_AUDIT_PLAYBOOK.md (audit type 2) in the
# workspace root.

skip_if_not_installed("move2")
skip_if_not_installed("sf")

library(move2)
library(sf)

read_f1 <- function(n = 300L) {
  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  fishers[mt_track_id(fishers) == "F1", ][seq_len(n), ]
}

# Same track in three CRSs: lon/lat (native), UTM 18N, and a local AEQD.
three_crs <- function(f) {
  list(
    longlat = f,
    utm     = sf::st_transform(f, sf::st_crs(32618)),
    aeqd    = sf::st_transform(f, move2::mt_aeqd_crs(f))
  )
}

flags_invariant <- function(fn) {
  reps <- three_crs(read_f1())
  out  <- lapply(reps, function(z)
    suppressWarnings(suppressMessages(fn(z))))
  identical(out$longlat, out$utm) && identical(out$longlat, out$aeqd)
}

# ---- outlier detectors: flags must be identical across input CRS ----------
# (the v0.4.1 CRS-stratification guarantee: geometry is canonicalised to a
# per-track AEQD regardless of input CRS, so flags are projection-independent.)

test_that("mt_clean_track is_outlier is invariant to input CRS", {
  expect_true(flags_invariant(function(z)
    mt_clean_track(z, remove = FALSE)$is_outlier))
})

test_that("mt_flag_outliers is_outlier is invariant to input CRS", {
  expect_true(flags_invariant(function(z) mt_flag_outliers(z)$is_outlier))
})

test_that("mt_flag_outliers_bridge is_outlier is invariant to input CRS", {
  expect_true(flags_invariant(function(z)
    mt_flag_outliers_bridge(z)$is_outlier))
})

test_that("mt_sequential_outliers is_outlier is invariant to input CRS", {
  expect_true(flags_invariant(function(z)
    mt_sequential_outliers(z)$is_outlier))
})

# ---- ud_outer_probability: query points auto-aligned to the UD's CRS ------
# Regression guard for the conceptual-integrity finding: a UD returned in one
# CRS queried with points in another silently returned NA / wrong cells.

test_that("ud_outer_probability is invariant to the query points' CRS", {
  skip_if_not_installed("terra")
  f   <- read_f1(200L)
  ud  <- suppressMessages(mt_dbbmm_ud(
    sf::st_transform(f, move2::mt_aeqd_crs(f)),
    location_error = 20, dim_size = 80, ext = 1.5))
  op_native <- ud_outer_probability(
    sf::st_transform(f, sf::st_crs(terra::crs(ud))), ud)  # already in UD CRS
  op_longlat <- ud_outer_probability(f, ud)               # lon/lat -> auto-reproject
  expect_equal(op_native, op_longlat)
  expect_true(sum(!is.na(op_longlat)) > 0)
})

# ---- mt_thin_distance: invariant to the CRS's linear unit -----------------
# Regression guard for the distance-unit leak (mt_distance() now units="m").
# NB: this asserts unit-invariance (metres vs US-feet of the SAME projection),
# NOT lon/lat-vs-projected invariance -- lon/lat uses geodesic distance while a
# projected CRS uses planar distance, so those legitimately differ. Making
# thin_distance fully geodesic/planar-invariant is a separate design choice.

test_that("mt_thin_distance is invariant to the CRS's linear unit", {
  f      <- read_f1()
  utm_m  <- sf::st_transform(f, sf::st_crs(32618))
  utm_ft <- sf::st_transform(
    f, sf::st_crs("+proj=utm +zone=18 +datum=WGS84 +units=us-ft +no_defs"))
  sel_m  <- mt_thin_distance(utm_m,  distance = 200)$thin_selected
  sel_ft <- mt_thin_distance(utm_ft, distance = 200)$thin_selected
  expect_identical(sel_m, sel_ft)
})

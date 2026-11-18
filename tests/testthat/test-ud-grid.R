## UD grid: arguments, the kernel-bounds probe, and grids large enough.
##
## The dBBMM and dBGB kernels (ported from move) stop when a bridge's kernel
## window leaves the grid.  Since the release after 0.4.5 an automatic grid is enlarged until
## every window fits, a grid that already fits is used unchanged, and a
## template is never enlarged.  See ?mt_dbbmm_ud, section "Grid".

skip_if_not_installed("move2")
skip_if_not_installed("terra")

ns <- asNamespace("move2utils")

## The worked example of the move2utils Application Note: Leroy (fisher M1),
## thinned to 30 min, projected to a track-local AEQD.  With v0.4.5 defaults
## both UD calls stopped with a grid error.
leroy_vars <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      fishers <- move2::mt_read(move2::mt_example())
      fishers <- fishers[!sf::st_is_empty(fishers), ]
      leroy <- move2::filter_track_data(fishers, .track_id = "M1")
      thin <- mt_thin_time(leroy,
                           interval = as.difftime(30, units = "mins"),
                           tolerance = as.difftime(2, units = "mins"),
                           remove = TRUE)
      la <- sf::st_transform(thin, move2::mt_aeqd_crs(thin))
      cache <<- list(
        track = la,
        iso = mt_dbbmm_variance(la, window_size = 31, margin = 11,
                                location_error = 20),
        dir = mt_dbgb_variance(la, window_size = 31, margin = 11,
                               location_error = 20))
    }
    cache
  }
})

kin_of <- function(v) {
  if (inherits(v, "mt_dbbmm_variance")) {
    ns$.dbbmm_kernel_inputs(v, NULL, NULL)
  } else {
    ns$.dbgb_kernel_inputs(v, NULL, NULL)
  }
}

pts_of <- function(v) {
  td <- v$track_data
  sf::st_as_sf(data.frame(x = td$x, y = td$y), coords = c("x", "y"),
               crs = sf::st_crs(td$crs))
}

test_that("the Note's worked example runs as printed (dBBMM defaults, dBGB ext 0.75)", {
  skip_on_cran()
  lv <- leroy_vars()
  ud_iso <- mt_dbbmm_ud(lv$iso, verbose = FALSE)
  ud_dir <- mt_dbgb_ud(lv$dir, ext_factor = 0.75, verbose = FALSE)
  for (ud in list(ud_iso, ud_dir)) {
    expect_s4_class(ud, "SpatRaster")
    expect_equal(sum(terra::values(ud), na.rm = TRUE), 1, tolerance = 1e-12)
    ## dim_size = 100 is kept on the longest side; cells are square
    expect_equal(max(dim(ud)[1:2]), 100)
    expect_equal(terra::res(ud)[1], terra::res(ud)[2], tolerance = 1e-9)
  }
  expect_true(ns$.ud_probe(kin_of(lv$iso), ud_iso)$safe)
  expect_true(ns$.ud_probe(kin_of(lv$dir), ud_dir)$safe)
})

test_that("enlarging is reported under verbose", {
  skip_on_cran()
  lv <- leroy_vars()
  expect_message(mt_dbbmm_ud(lv$iso, verbose = TRUE),
                 class = "move2utils_ud_grid_enlarged")
})

test_that("a grid that is already large enough is used unchanged", {
  skip_on_cran()
  lv <- leroy_vars()
  for (v in list(lv$iso, lv$dir)) {
    ga <- ns$.ud_grid_args(NULL, NULL, 100, 1.5, NULL, TRUE)
    old <- ns$.make_raster(pts_of(v), dim_size = 100, ext = 1.5)
    expect_true(ns$.ud_probe(kin_of(v), old)$safe)
    new <- ns$.ud_auto_grid(pts_of(v), ga, list(kin_of(v)))
    expect_identical(as.vector(terra::ext(new)), as.vector(terra::ext(old)))
    expect_identical(dim(new), dim(old))
  }
})

test_that("res keeps the cell size and dim_size keeps the cell count", {
  skip_on_cran()
  lv <- leroy_vars()
  for (v in list(lv$iso, lv$dir)) {
    ## both too small at ext_factor = 0.1
    ga_res <- ns$.ud_grid_args(NULL, 50, NULL, 0.1, NULL, TRUE)
    g_res <- ns$.ud_auto_grid(pts_of(v), ga_res, list(kin_of(v)))
    expect_equal(terra::res(g_res), c(50, 50), tolerance = 1e-9)
    expect_true(ns$.ud_probe(kin_of(v), g_res)$safe)

    ga_dim <- ns$.ud_grid_args(NULL, NULL, 80, 0.1, NULL, TRUE)
    g_dim <- ns$.ud_auto_grid(pts_of(v), ga_dim, list(kin_of(v)))
    expect_equal(max(dim(g_dim)[1:2]), 80)
    expect_true(ns$.ud_probe(kin_of(v), g_dim)$safe)
  }
})

test_that("a grid from res or dim_size is capped even when not enlarged", {
  skip_on_cran()
  lv <- leroy_vars()
  old <- options(move2utils.ud_max_cells = NULL)
  on.exit(options(old), add = TRUE)
  v <- lv$iso
  for (ga in list(ns$.ud_grid_args(NULL, 100, NULL, 1.5, NULL, TRUE),
                  ns$.ud_grid_args(NULL, NULL, 200, 1.5, NULL, TRUE))) {
    g <- ns$.ud_auto_grid(pts_of(v), ga, list(kin_of(v)))
    ## built as is (large enough without enlargement)
    expect_identical(dim(g)[1:2],
                     dim(if (is.null(ga$res)) {
                       ns$.make_raster(pts_of(v), dim_size = ga$dim_size,
                                       ext = ga$ext_factor)
                     } else {
                       ns$.make_raster(pts_of(v), cell_size = ga$res,
                                       ext = ga$ext_factor)
                     })[1:2])
    options(move2utils.ud_max_cells = terra::ncell(g) - 1)
    expect_error(ns$.ud_auto_grid(pts_of(v), ga, list(kin_of(v))),
                 class = "move2utils_ud_grid_too_large")
    options(move2utils.ud_max_cells = NULL)
  }
})

test_that("an enlarged res grid is capped by move2utils.ud_max_cells, for both models", {
  skip_on_cran()
  lv <- leroy_vars()
  old <- options(move2utils.ud_max_cells = NULL)
  on.exit(options(old), add = TRUE)
  for (v in list(lv$iso, lv$dir)) {
    ga <- ns$.ud_grid_args(NULL, 50, NULL, 0.1, NULL, TRUE)
    g <- ns$.ud_auto_grid(pts_of(v), ga, list(kin_of(v)))
    n <- terra::ncell(g)
    ## just under the size of the grid that is needed: stop before allocating
    options(move2utils.ud_max_cells = n - 1)
    expect_error(ns$.ud_auto_grid(pts_of(v), ga, list(kin_of(v))),
                 class = "move2utils_ud_grid_too_large")
    ## at exactly that size: the same grid as without a cap
    options(move2utils.ud_max_cells = n)
    expect_identical(dim(ns$.ud_auto_grid(pts_of(v), ga, list(kin_of(v)))),
                     dim(g))
    options(move2utils.ud_max_cells = NULL)
  }
  options(move2utils.ud_max_cells = "many")
  expect_error(ns$.ud_max_cells(), class = "move2utils_ud_bad_option")
})

test_that("the kernel runs cleanly on every grid the probe calls safe", {
  ## The raw kernel is only ever called on grids the probe has cleared: on an
  ## unchecked grid it can write outside the grid before its own check stops
  ## it, which is the defect the probe exists to prevent.
  skip_on_cran()
  lv <- leroy_vars()
  kin <- kin_of(lv$iso)
  for (e in c(1.5, 2)) {
    g <- ns$.make_raster(pts_of(lv$iso), dim_size = 60, ext = e)
    expect_true(ns$.ud_probe(kin, g)$safe, info = sprintf("ext = %g", e))
    gv <- ns$.ud_grid_vectors(g, "dbbmm")
    raw <- .Call("dbbmm2_omp", kin$x, kin$y, kin$s1, kin$t, kin$le, gv$x,
                 gv$y, kin$time_step, kin$ext, kin$interest,
                 PACKAGE = "move2utils")
    ud <- mt_dbbmm_ud(lv$iso, dim_size = 60, ext_factor = e, verbose = FALSE)
    ## lay the kernel's matrix onto the grid exactly as mt_dbbmm_ud() does
    terra::values(g) <- raw / sum(raw)
    expect_identical(terra::values(g), terra::values(ud))
  }
  ## and a grid the probe rejects never reaches a kernel through the API
  small <- ns$.make_raster(pts_of(lv$iso), dim_size = 60, ext = 0.2)
  expect_false(ns$.ud_probe(kin, small)$safe)
  expect_error(mt_dbbmm_ud(lv$iso, raster = small, verbose = FALSE),
               class = "move2utils_ud_raster_too_small")
})

test_that("a template that is too small stops, and the suggested extend() fixes it", {
  skip_on_cran()
  lv <- leroy_vars()
  small <- ns$.make_raster(pts_of(lv$iso), cell_size = 100, ext = 0.1)
  err <- tryCatch(mt_dbbmm_ud(lv$iso, raster = small, verbose = FALSE),
                  error = function(e) e)
  expect_s3_class(err, "move2utils_ud_raster_too_small")
  call_txt <- regmatches(conditionMessage(err),
                         regexpr("terra::extend\\(raster, [^\n]*\\)",
                                 conditionMessage(err)))
  expect_length(call_txt, 1L)
  raster <- small
  bigger <- eval(parse(text = call_txt))
  expect_equal(terra::res(bigger), terra::res(small))
  ud <- mt_dbbmm_ud(lv$iso, raster = bigger, verbose = FALSE)
  expect_equal(sum(terra::values(ud), na.rm = TRUE), 1, tolerance = 1e-12)
  ## the template is used as given, never enlarged
  expect_identical(as.vector(terra::ext(ud)), as.vector(terra::ext(bigger)))
})

test_that("grid arguments are validated", {
  skip_on_cran()
  lv <- leroy_vars()
  tmpl <- ns$.make_raster(pts_of(lv$iso), cell_size = 100, ext = 1.5)
  expect_error(mt_dbbmm_ud(lv$iso, raster = tmpl, res = 100),
               class = "move2utils_ud_grid_overspecified")
  expect_error(mt_dbbmm_ud(lv$iso, res = 100, dim_size = 50),
               class = "move2utils_ud_grid_overspecified")
  expect_error(mt_dbgb_ud(lv$dir, raster = tmpl, dim_size = 50),
               class = "move2utils_ud_grid_overspecified")
  expect_error(mt_dbbmm_ud(lv$iso, res = c(100, 50)),
               class = "move2utils_ud_bad_res")
  expect_error(mt_dbbmm_ud(lv$iso, dim_size = 1),
               class = "move2utils_ud_bad_dim_size")
  expect_error(mt_dbbmm_ud(lv$iso, dim_size = 50.5),
               class = "move2utils_ud_bad_dim_size")
  expect_error(mt_dbbmm_ud(lv$iso, ext_factor = -1),
               class = "move2utils_ud_bad_ext_factor")
  expect_error(mt_dbbmm_ud(lv$iso, ext = terra::ext(tmpl)),
               class = "move2utils_ud_ext_is_extent")
  expect_error(mt_dbgb_ud(lv$dir, ext_factor = terra::ext(tmpl)),
               class = "move2utils_ud_ext_is_extent")
  expect_error(mt_dbbmm_ud(lv$iso, ext = 1, ext_factor = 1),
               class = "move2utils_ud_ext_twice")
  expect_error(mt_dbbmm_ud(lv$iso, raster = "100"),
               class = "move2utils_ud_bad_raster")
})

test_that("the deprecated forms warn and give identical results", {
  skip_on_cran()
  lv <- leroy_vars()
  expect_warning(a <- mt_dbbmm_ud(lv$iso, ext = 1.5, dim_size = 50,
                                         verbose = FALSE),
                 class = "move2utils_deprecated_ext")
  b <- mt_dbbmm_ud(lv$iso, ext_factor = 1.5, dim_size = 50, verbose = FALSE)
  expect_identical(terra::values(a), terra::values(b))
  expect_identical(as.vector(terra::ext(a)), as.vector(terra::ext(b)))

  expect_warning(c1 <- mt_dbgb_ud(lv$dir, raster = 200, ext_factor = 1.5,
                                  verbose = FALSE),
                 class = "move2utils_deprecated_raster_numeric")
  d1 <- mt_dbgb_ud(lv$dir, res = 200, ext_factor = 1.5, verbose = FALSE)
  expect_identical(terra::values(c1), terra::values(d1))
})

test_that("the multi-track common grid holds every track's kernel windows", {
  skip_on_cran()
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  two <- fishers[move2::mt_track_id(fishers) %in% c("M1", "M4"), ]
  two <- mt_thin_time(two, interval = as.difftime(60, units = "mins"),
                      tolerance = as.difftime(5, units = "mins"),
                      remove = TRUE)
  two <- sf::st_transform(two, move2::mt_aeqd_crs(two))
  vars <- mt_dbbmm_variance(two, window_size = 31, margin = 11,
                            location_error = 20)
  stk <- mt_dbbmm_ud(vars, ext_factor = 0.1, dim_size = 60, verbose = FALSE)
  expect_equal(terra::nlyr(stk), 2L)
  for (nm in names(vars)) {
    expect_true(ns$.ud_probe(kin_of(vars[[nm]]), stk)$safe, info = nm)
    expect_equal(sum(terra::values(stk[[nm]]), na.rm = TRUE), 1,
                 tolerance = 1e-12)
  }
})

test_that("degenerate grids and inputs stop before a kernel runs", {
  skip_on_cran()
  lv <- leroy_vars()
  ## a one-row template: the kernels would read past the grid coordinates
  e <- as.vector(terra::ext(ns$.make_raster(pts_of(lv$iso), cell_size = 100,
                                             ext = 1.5)))
  one_row <- terra::rast(nrows = 1, ncols = 50, xmin = e[1], xmax = e[2],
                         ymin = e[3], ymax = e[4],
                         crs = sf::st_crs(lv$iso$track_data$crs)$wkt)
  expect_false(ns$.ud_probe(kin_of(lv$iso), one_row)$safe)
  expect_error(mt_dbbmm_ud(lv$iso, raster = one_row, verbose = FALSE),
               class = "move2utils_ud_raster_too_small")
  expect_error(mt_dbgb_ud(lv$dir, raster = one_row, verbose = FALSE),
               class = "move2utils_ud_raster_too_small")
  ## a time step that is not positive would never end
  expect_error(mt_dbbmm_ud(lv$iso, time_step = 0, verbose = FALSE),
               class = "move2utils_ud_bad_time_step")
  expect_error(mt_dbgb_ud(lv$dir, time_step = -1, verbose = FALSE),
               class = "move2utils_ud_bad_time_step")
  ## dBGB with no bridge of interest
  one <- lv$dir
  one$seg_interest[] <- FALSE
  expect_error(mt_dbgb_ud(one, verbose = FALSE),
               class = "move2utils_ud_too_few_fixes")
  ## dim_size 2 cannot hold the kernel margin; 3 and up can
  expect_error(mt_dbbmm_ud(lv$iso, dim_size = 2, ext_factor = 0, verbose = FALSE),
               class = "move2utils_ud_dim_size_too_small")
  for (d in 3:5) {
    ga <- ns$.ud_grid_args(NULL, NULL, d, 0, NULL, TRUE)
    g <- ns$.ud_auto_grid(pts_of(lv$iso), ga, list(kin_of(lv$iso)))
    expect_equal(max(dim(g)[1:2]), d, info = sprintf("dim_size = %d", d))
    expect_true(ns$.ud_probe(kin_of(lv$iso), g)$safe)
  }
})

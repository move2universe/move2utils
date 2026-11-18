#' Internal UD utilities
#'
#' Internal imports and dynamic library registration for the dBBMM/dBGB
#' variance and UD machinery.
#'
#' @importFrom sf st_is_longlat st_bbox st_as_sf st_crs
#' @importFrom move2 mt_n_tracks mt_is_move2 mt_time mt_track_id mt_time_lags
#' @importFrom terra xFromCol yFromRow ncol nrow ncell
#' @importFrom units set_units as_units
#' @importFrom stats aggregate optim
#' @useDynLib move2utils, .registration = TRUE
#' @name move2utils-ud-internals
#' @keywords internal
NULL

## Internal: factory for the inline `say` closure used by chatty
## primitives (mt_clean_track, mt_flag_outliers, mt_flag_outliers_dbgb).
## Avoids open-coding `say <- function(...) if (!silent) message(...)`
## at every call site.
## @keywords internal
.say <- function(silent) {
  function(...) if (!silent) message(...)
}

#' Extract coordinates, timestamps, and time lags from a move2 object
#'
#' Internal helper that validates and extracts the numeric vectors
#' needed by the dBBMM/dBGB algorithms. Performs all input validation.
#'
#' @param x A `move2` object (single track).
#' @return A list with components `x`, `y`, `time_mins`, `n_locs`.
#' @keywords internal
.extract_track_data <- function(x) {
  .validate_move2(x)

  if (mt_n_tracks(x) > 1) {
    rlang::abort(
      "Internal error: .extract_track_data called with multiple tracks.",
      class = "move2utils_internal_extract_track_data_multi_track")
  }

  coords <- st_coordinates(x)
  if (any(!is.finite(coords))) {
    rlang::abort("Input contains non-finite coordinates (NA, NaN, or Inf).",
                 class = "move2utils_input_nonfinite_coords")
  }

  ts_numeric <- as.numeric(mt_time(x))
  if (is.unsorted(ts_numeric, strictly = FALSE)) {
    rlang::abort("Timestamps are not in chronological order.",
                 class = "move2utils_input_unsorted_timestamps")
  }
  if (any(diff(ts_numeric) == 0)) {
    rlang::abort(paste0(
      "Input contains duplicate timestamps. Remove duplicates first: ",
      "move2::mt_filter_unique(x)"),
      class = "move2utils_input_duplicate_timestamps")
  }

  time_mins <- ts_numeric / 60

  list(
    x = coords[, 1],
    y = coords[, 2],
    time_mins = time_mins,
    n_locs = nrow(coords)
  )
}

#' Expand location error to match number of locations
#' @keywords internal
.expand_loc_error <- function(location_error, n_locs) {
  if (length(location_error) == 1) {
    location_error <- rep(location_error, n_locs)
  }
  if (length(location_error) != n_locs) {
    rlang::abort(
      "location_error must be length 1 or equal to the number of locations.",
      class = "move2utils_location_error_bad_vector_length")
  }
  if (any(is.na(location_error))) {
    rlang::abort("location_error must not contain NAs.",
                 class = "move2utils_location_error_has_na")
  }
  if (any(location_error < 0)) {
    rlang::abort("location_error must be non-negative.",
                 class = "move2utils_location_error_negative")
  }
  location_error
}

#' Validate a move2 object (common checks that apply regardless of track count)
#' @keywords internal
.validate_move2 <- function(x) {
  if (!mt_is_move2(x)) {
    rlang::abort("Input must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  if (st_is_longlat(x)) {
    rlang::abort(paste0(
      "Cannot use longitude/latitude coordinates. ",
      "Transform to a projected CRS first, e.g.: ",
      "sf::st_transform(x, move2::mt_aeqd_crs(x))"),
      class = "move2utils_input_longlat_required_projected")
  }
  empties <- st_is_empty(x)
  if (any(empties)) {
    rlang::abort(
      sprintf("Input contains %d empty geometries (failed GPS fixes). Remove them first: x <- x[!sf::st_is_empty(x), ]",
              sum(empties)),
      class = "move2utils_input_has_empty_geometries")
  }
  invisible(TRUE)
}

#' Ensure a move2 object is in a metric (projected) CRS, reprojecting from
#' longitude/latitude to a local AEQD if necessary.  The benefit of taking a
#' move2 object is that the CRS is known, so functions that need metric units
#' (variance, UD) can reproject internally rather than forcing the user to do
#' it.  Non-move2 input is returned unchanged for the caller's own validation
#' to reject.  Emits an informative message when it reprojects (suppressible by
#' the caller wrapping in suppressMessages()).
#' @keywords internal
.ensure_projected <- function(x, what = "this computation") {
  if (mt_is_move2(x)) {
    crs <- sf::st_crs(x)
    if (!is.na(crs) && sf::st_is_longlat(x)) {
      message(sprintf(
        paste0("Input is in longitude/latitude.  Auto-projecting to a local ",
               "AEQD CRS for %s; distances and variances are in metres."),
        what))
      x <- sf::st_transform(x, move2::mt_aeqd_crs(x))
    }
  }
  x
}

#' Reject empty geometries at a function's entry.  There is no outlier to
#' identify in a fix with no location, so empty geometry must be removed
#' upstream rather than silently dropped.  Errors with the shared
#' `move2utils_input_has_empty_geometries` class.
#' @keywords internal
.reject_empty_geometry <- function(x, what = "this computation") {
  if (!mt_is_move2(x)) return(invisible(x))
  empty <- sf::st_is_empty(x)
  if (any(empty)) {
    rlang::abort(
      sprintf(paste0(
        "%d fix(es) have empty geometry. Remove them before %s, e.g. ",
        "`x <- x[!sf::st_is_empty(x), ]` or `mt_filter_gps_quality(x)`. ",
        "There is no outlier to identify in a fix with no location."),
        sum(empty), what),
      class = "move2utils_input_has_empty_geometries")
  }
  invisible(x)
}

#' Is a CRS suitable for running the Brownian-bridge (dBBMM/dBGB) math
#' directly?  The kernel needs Euclidean distances in metres, so a CRS is
#' suitable iff it is projected (not longitude/latitude) AND its linear unit
#' is the metre.  We deliberately do NOT second-guess a metric projection the
#' user chose (e.g. flagging high-latitude web-Mercator distortion): trusting
#' the user's projection is the contract for the legacy UD layer.  A non-metre
#' linear unit (e.g. US-feet) returns FALSE and routes through the AEQD path --
#' correct numerics, just a resample at the end.
#' @keywords internal
.crs_is_suitable <- function(crs) {
  crs <- tryCatch(sf::st_crs(crs), error = function(e) sf::NA_crs_)
  if (is.na(crs)) return(FALSE)
  if (isTRUE(sf::st_is_longlat(crs))) return(FALSE)
  u <- crs$units_gdal
  if (is.null(u) || is.na(u)) return(FALSE)
  grepl("met", u, ignore.case = TRUE)   # "metre" / "meter"
}

#' Warp a UD raster (cell value = probability mass, summing to 1) onto a
#' target grid and renormalise.  `target` is either a `terra::SpatRaster`
#' template (the result lands on its exact grid -- used to align a UD to a
#' user-supplied environmental layer) or a CRS (terra picks a grid).
#' `terra::project()` resamples, so the warped mass no longer sums to 1 and
#' edge cells become NA / can overshoot slightly negative; both are corrected
#' here, per layer, so the returned UD again honours the "mass sums to 1"
#' invariant that downstream consumers (ud_outer_probability, manual cell
#' reads) rely on.  (`ud_volume()` renormalises independently, so the contour
#' representation is safe regardless; this protects the raw UD raster.)
#' @keywords internal
.warp_ud_to_target <- function(r, target) {
  tgt <- if (inherits(target, "SpatRaster")) target else sf::st_crs(target)$wkt
  out <- terra::project(r, tgt)
  vals <- terra::values(out)
  vals[!is.finite(vals) | vals < 0] <- 0
  for (j in seq_len(ncol(vals))) {
    s <- sum(vals[, j], na.rm = TRUE)
    if (is.finite(s) && s > 0) vals[, j] <- vals[, j] / s
  }
  terra::values(out) <- vals
  out
}

#' Resolve, for a single-track UD computation, (1) the raster grid the kernel
#' runs on -- always in the variance object's metric *compute* CRS
#' (`td$crs`) -- and (2) the post-compute warp target, i.e. where the UD is
#' returned.  Implements the CRS-stratification contract for the legacy UD
#' layer (see `mt_dbbmm_ud` / `mt_dbgb_ud`): the *target* CRS is the
#' user-supplied template's CRS if given, else `target_crs` (passed down from
#' the `.move2` entry), else the variance object's recorded `orig_crs`.  When
#' the target equals the compute CRS the UD is computed directly on it (no
#' warp, lossless -- the "suitable CRS" fast path); otherwise the result is
#' warped onto the target at the end.  An explicit template in another CRS is
#' reprojected into the compute CRS for the kernel and the result is warped
#' back onto the template's exact grid.  Returns `list(grid, warp_to)` where
#' `warp_to` is a SpatRaster template, a CRS, or NULL (no warp needed).
#' A template is checked against the kernel windows of `kins` and never
#' enlarged; an automatic grid is built large enough (`.ud_auto_grid()`).
#' @keywords internal
.ud_resolve_grid <- function(grid_args, kins, td, target_crs, restore,
                             verbose = FALSE) {
  raster <- grid_args$raster
  compute_crs <- sf::st_crs(td$crs)
  if (!is.null(target_crs)) {
    target <- sf::st_crs(target_crs)
  } else if (inherits(raster, "SpatRaster")) {
    target <- sf::st_crs(terra::crs(raster))
  } else if (!is.null(td$orig_crs)) {
    target <- sf::st_crs(td$orig_crs)
  } else {
    target <- compute_crs
  }

  if (inherits(raster, "SpatRaster")) {
    if (sf::st_crs(terra::crs(raster)) != compute_crs) {
      ## explicit environmental grid in another CRS: compute on its
      ## reprojection into the metric compute CRS, warp the result back onto
      ## the exact env grid the user handed in.
      grid <- terra::project(raster, compute_crs$wkt)
      .ud_check_template(grid, kins, reprojected = TRUE)
      return(list(grid = grid, warp_to = raster))
    }
    .ud_check_template(raster, kins, reprojected = FALSE)
    warp_to <- if (isTRUE(restore) && target != compute_crs) target else NULL
    return(list(grid = raster, warp_to = warp_to))
  }

  ## auto grid (from `res` or `dim_size`) built in the metric compute CRS
  pts <- sf::st_as_sf(data.frame(x = td$x, y = td$y),
                      coords = c("x", "y"), crs = compute_crs)
  grid <- .ud_auto_grid(pts, grid_args, kins, verbose = verbose)
  warp_to <- if (isTRUE(restore) && target != compute_crs) target else NULL
  list(grid = grid, warp_to = warp_to)
}

## ---------------------------------------------------------------------------
## UD grid: arguments, the kernels' bounds, and a grid large enough.
##
## dbbmm2_omp() and bgb_omp() spread each bridge time step over a window of 4
## (dBBMM) or 5 (dBGB) standard deviations and stop when that window leaves
## the grid.  `ud_grid_probe` (src/ud_grid_probe.c) walks the same time steps
## with the same arithmetic and reports the window envelope and the index
## extremes on a candidate grid, so the grid can be checked, and if need be
## enlarged, before a kernel runs.  The kernels are unchanged.
## ---------------------------------------------------------------------------

#' Normalise the grid arguments of mt_dbbmm_ud() / mt_dbgb_ud()
#'
#' Exactly one of `raster` (a template SpatRaster), `res` and `dim_size` may
#' be given; with none, `dim_size = 100`.  The deprecated forms (`ext`, and a
#' number passed as `raster`) are mapped onto `ext_factor` and `res` with a
#' warning.
#' @keywords internal
.ud_grid_args <- function(raster, res, dim_size, ext_factor, ext,
                          ext_factor_given) {
  if (inherits(ext, "SpatExtent") || inherits(ext_factor, "SpatExtent")) {
    rlang::abort(paste0(
      "`ext_factor` (formerly `ext`) is a padding factor, not an extent: the ",
      "track's bounding box is widened by this fraction of its width and ",
      "height on each side. To compute the UD over a given extent, pass a ",
      "template raster instead, e.g. ",
      "`raster = terra::rast(<extent>, resolution = <cell size>, crs = <crs>)`."),
      class = "move2utils_ud_ext_is_extent")
  }
  if (!is.null(ext)) {
    if (isTRUE(ext_factor_given)) {
      rlang::abort("Give `ext_factor` only. `ext` is its deprecated former name.",
                   class = "move2utils_ud_ext_twice")
    }
    rlang::warn(paste0(
      "`ext` is deprecated; use `ext_factor` (same meaning: the fraction by ",
      "which the track's bounding box is widened). Renamed because ",
      "`terra::ext()` denotes an extent."),
      class = c("move2utils_deprecated_ext", "move2utils_deprecated"))
    ext_factor <- ext
  }
  if (!is.numeric(ext_factor) || length(ext_factor) != 1L ||
      !is.finite(ext_factor) || ext_factor < 0) {
    rlang::abort("`ext_factor` must be a single non-negative number.",
                 class = "move2utils_ud_bad_ext_factor")
  }
  if (is.numeric(raster)) {
    if (!is.null(res)) {
      rlang::abort(paste0(
        "Give the cell size as `res` only. A number passed as `raster` is ",
        "the deprecated way of giving the same thing."),
        class = "move2utils_ud_grid_overspecified")
    }
    rlang::warn(paste0(
      "Passing a cell size as `raster` is deprecated; use `res`. `raster` ",
      "now takes a template SpatRaster only."),
      class = c("move2utils_deprecated_raster_numeric",
                "move2utils_deprecated"))
    res <- raster
    raster <- NULL
  }
  if (!is.null(raster) && !inherits(raster, "SpatRaster")) {
    rlang::abort("`raster` must be a terra SpatRaster (a template grid) or NULL.",
                 class = "move2utils_ud_bad_raster")
  }
  given <- c(raster = !is.null(raster), res = !is.null(res),
             dim_size = !is.null(dim_size))
  if (sum(given) > 1L) {
    rlang::abort(paste0(
      "Give only one of `raster` (a template), `res` (a cell size) or ",
      "`dim_size` (cells along the longest side); got ",
      paste0("`", names(given)[given], "`", collapse = " and "), "."),
      class = "move2utils_ud_grid_overspecified")
  }
  if (!is.null(res) && (!is.numeric(res) || length(res) != 1L ||
                        !is.finite(res) || res <= 0)) {
    rlang::abort(paste0(
      "`res` must be a single positive number: the cell size in map units ",
      "of the compute CRS. Cells are square."),
      class = "move2utils_ud_bad_res")
  }
  if (!is.null(dim_size) && (!is.numeric(dim_size) || length(dim_size) != 1L ||
                             !is.finite(dim_size) || dim_size < 2 ||
                             dim_size != round(dim_size))) {
    rlang::abort(paste0(
      "`dim_size` must be a single whole number of at least 2: the number of ",
      "cells along the longest side. Cells are square."),
      class = "move2utils_ud_bad_dim_size")
  }
  if (sum(given) == 0L) dim_size <- 100
  list(raster = raster, res = res, dim_size = dim_size,
       ext_factor = ext_factor)
}

#' Resolve a location-error input against a variance object
#' @keywords internal
.ud_location_error <- function(td, location_error) {
  if (is.null(location_error)) {
    ## Pull from the stored per-fix vector if present; fall back to 0 for
    ## variance objects fit before the storage field existed (backward
    ## compatibility).
    if (!is.null(td$location_error)) td$location_error else rep(0, td$n_locs)
  } else {
    .expand_loc_error(location_error, td$n_locs)
  }
}

#' Stop before a kernel is called on inputs it cannot handle
#'
#' Both kernels read `t[k+1]` while walking segments, so they need at least
#' two fixes (of interest, for dBGB), and they step time by `time_step`, so
#' a step that is not positive never ends.
#' @keywords internal
.ud_check_kernel_inputs <- function(kin) {
  if (length(kin$x) < 2L) {
    rlang::abort(paste0(
      "A UD needs at least two fixes",
      if (kin$model == "dbgb") " of interest (see `seg_interest`)" else "",
      "; got ", length(kin$x), "."),
      class = "move2utils_ud_too_few_fixes")
  }
  if (kin$model == "dbgb" && !any(kin$interest[-length(kin$interest)])) {
    rlang::abort(paste0(
      "A UD needs at least one bridge of interest; every bridge is ",
      "excluded (see `seg_interest` and mt_mask_segments())."),
      class = "move2utils_ud_too_few_fixes")
  }
  if (!is.numeric(kin$time_step) || length(kin$time_step) != 1L ||
      !is.finite(kin$time_step) || kin$time_step <= 0) {
    rlang::abort(if (isTRUE(kin$time_step_derived)) {
      paste0("The default `time_step` is a fraction of the shortest time ",
             "lag, which is not positive here (duplicate or unsorted ",
             "timestamps?). Give `time_step` in minutes.")
    } else {
      "`time_step` must be a single positive number (minutes)."
    }, class = "move2utils_ud_bad_time_step")
  }
  kin
}

#' The exact inputs of dbbmm2_omp(), built once for the probe and the kernel
#' @keywords internal
.dbbmm_kernel_inputs <- function(object, location_error, time_step) {
  td <- object$track_data
  location_error <- .ud_location_error(td, location_error)
  time_lag <- c(diff(td$time_mins), 0)
  time_step_derived <- is.null(time_step)
  if (time_step_derived) {
    time_step <- min(time_lag[-length(time_lag)]) / 15
  }
  ## NB: dbbmm2_omp's third argument is the per-fix BM variance (sigma^2),
  ## not the standard deviation -- consistent with the legacy
  ## `move::brownian.motion.variance.dyn` -> @means slot (also variance
  ## scale) and with mt_motion_variance() which returns variance for both
  ## dBBMM and dBGB (the dBGB accessor squares the internally-stored sigmas).
  variance <- c(object$variance, 0)
  variance[is.na(variance)] <- 0
  .ud_check_kernel_inputs(list(
    model = "dbbmm", x = td$x, y = td$y, s1 = variance, s2 = variance,
    t = td$time_mins - min(td$time_mins), le = location_error,
    time_step = time_step, ext = 4, interest = object$interest,
    time_lag = time_lag, time_step_derived = time_step_derived))
}

#' The exact inputs of bgb_omp(), built once for the probe and the kernel
#' @keywords internal
.dbgb_kernel_inputs <- function(object, location_error, time_step) {
  td <- object$track_data
  n <- td$n_locs
  location_error <- .ud_location_error(td, location_error)
  # Two masks, kept apart.  The core is the run of fixes whose variance
  # window was fully covered (n_estim at its maximum); only these enter the
  # kernel, as in move.  It is derived from n_estim, never from
  # seg_interest, so that masking can never drop a fix.  seg_interest is
  # the per-bridge include mask: seg_interest[i] = TRUE means "include the
  # bridge from location i to location i+1 in the UD", and the kernel skips
  # bridges marked FALSE exactly as dbbmm2_omp() does with `interest` (see
  # mt_mask_segments()).  Unmasked, every core bridge counts and the result
  # is the one move gives.  (Before 0.4.6 the kernel was subset to
  # seg_interest, so masking dropped the fix and bridged across it.)
  ne <- object$n_estim
  core <- if (is.null(ne) || all(is.na(ne))) {
    object$seg_interest
  } else {
    !is.na(ne) & ne == max(ne, na.rm = TRUE)
  }
  keep_bridge <- object$seg_interest[core]
  t_mins <- td$time_mins
  time_step_derived <- is.null(time_step)
  if (time_step_derived) {
    time_step <- min(diff(t_mins)) / 20.1
  }
  para_sd <- object$para_sd
  orth_sd <- object$orth_sd
  para_sd[is.na(para_sd)] <- 0
  orth_sd[is.na(orth_sd)] <- 0
  .ud_check_kernel_inputs(list(
    model = "dbgb", x = td$x[core], y = td$y[core],
    s1 = para_sd[core], s2 = orth_sd[core], t = t_mins[core],
    le = rep(location_error, length.out = n)[core],
    time_step = time_step, ext = 5, interest = keep_bridge,
    time_step_derived = time_step_derived))
}

#' Grid-coordinate vectors exactly as each kernel receives them
#' @keywords internal
.ud_grid_vectors <- function(grid, model) {
  x <- terra::xFromCol(grid, 1:terra::ncol(grid))
  y <- if (model == "dbbmm") {
    terra::yFromRow(grid, terra::nrow(grid):1)
  } else {
    sort(unique(terra::yFromRow(grid, 1:terra::nrow(grid))))
  }
  list(x = x, y = y)
}

#' Run the bounds probe for one track on one grid
#'
#' `safe` is TRUE when the kernel's window never leaves the grid, including
#' the equality cases the kernels' own check lets through (`xEnd == nX`,
#' `yEnd == nY`) at which they would index one past the grid.  A grid needs
#' at least two rows and two columns: the kernels take the cell size from
#' the first two grid coordinates.
#' @keywords internal
.ud_probe <- function(kin, grid) {
  gv <- .ud_grid_vectors(grid, kin$model)
  nx <- length(gv$x)
  ny <- length(gv$y)
  if (nx < 2L || ny < 2L) {
    ## the probe cannot run on such a grid, and neither can the kernels;
    ## report it as too small on the short side(s), with no envelope
    return(list(env = c(Inf, -Inf, Inf, -Inf),
                sides = c(left = nx < 2L, right = nx < 2L,
                          bottom = ny < 2L, top = ny < 2L),
                safe = FALSE))
  }
  r <- .Call("ud_grid_probe",
             if (kin$model == "dbbmm") 0L else 1L,
             kin$x, kin$y, kin$s1, kin$s2, kin$t, kin$le,
             gv$x, gv$y, kin$time_step, kin$ext, kin$interest)
  idx <- r[5:8]
  list(env = r[1:4],
       sides = c(left   = idx[1] < 0,
                 right  = idx[2] > nx - 1,
                 bottom = idx[3] > ny - 1,
                 top    = idx[4] < 0),
       safe = idx[1] >= 0 && idx[2] <= nx - 1 && idx[3] <= ny - 1 &&
              idx[4] >= 0)
}

#' Union of the window envelopes of several probes
#' @keywords internal
.ud_envelope_union <- function(probes) {
  e <- vapply(probes, function(p) p$env, numeric(4))
  c(min(e[1, ]), max(e[2, ]), min(e[3, ]), max(e[4, ]))
}

#' The extent a grid of square cells of size `cell` needs so that every
#' kernel window stays inside it, given the window envelope `env`.
#'
#' The kernels reach cells whose centres lie within the window, so the
#' envelope has to lie between the first and last cell centres: half a cell
#' beyond it left, right and top.  The y window is one row lower than the x
#' window (rows floor(v)-1 .. ceil(v')-1, inherited from move), so the
#' bottom needs one and a half cells.  `extra` widens every side by further
#' whole cells (used only if floating-point rounding defeats the minimum).
#' @keywords internal
.ud_needed_range <- function(box, env, cell, extra = 0) {
  c(xmin = min(box[["xmin"]], env[1] - (0.5 + extra) * cell),
    xmax = max(box[["xmax"]], env[2] + (0.5 + extra) * cell),
    ymin = min(box[["ymin"]], env[3] - (1.5 + extra) * cell),
    ymax = max(box[["ymax"]], env[4] + (0.5 + extra) * cell))
}

#' A raster of square cells covering `range`, snapped outward symmetrically
#'
#' The same construction as `.make_raster()`, except that the cell counts
#' are taken with a relative tolerance, so that a side that is a whole number
#' of cells long does not gain a spurious extra cell from rounding.
#' @keywords internal
.raster_from_range <- function(range, cell, crs_wkt) {
  x_range <- range[["xmax"]] - range[["xmin"]]
  y_range <- range[["ymax"]] - range[["ymin"]]
  nc <- max(1, ceiling(x_range / cell - 1e-9))
  nr <- max(1, ceiling(y_range / cell - 1e-9))
  xpad <- (nc * cell - x_range) / 2
  ypad <- (nr * cell - y_range) / 2
  terra::rast(nrows = nr, ncols = nc,
              xmin = range[["xmin"]] - xpad, xmax = range[["xmax"]] + xpad,
              ymin = range[["ymin"]] - ypad, ymax = range[["ymax"]] + ypad,
              crs = crs_wkt)
}

#' An automatic grid that is large enough for every kernel window
#'
#' The grid `ext_factor` gives (exactly as in v0.4.5 and `move`) is used
#' unchanged whenever it is large enough, so those results are identical.
#' Otherwise the extent is that grid's padded bounding box combined with the
#' window envelope of every track in `kins`.  With `res`, the cell size is
#' kept and cells are added; with `dim_size`, the longest side keeps that many
#' cells and the cell size grows (iterated to a fixed point, since a larger
#' cell needs a wider margin).
#' @keywords internal
.ud_auto_grid <- function(pts, grid_args, kins, verbose = FALSE) {
  res <- grid_args$res
  dim_size <- grid_args$dim_size
  ext_factor <- grid_args$ext_factor
  grid <- if (!is.null(res)) {
    .make_raster(pts, cell_size = res, ext = ext_factor)
  } else {
    .make_raster(pts, dim_size = dim_size, ext = ext_factor)
  }
  ## every automatic grid is capped, not only an enlarged one: a terra
  ## raster holds no values yet, so this stops before anything is allocated
  .ud_check_grid_cells(grid, res = res, dim_size = dim_size)
  probes <- lapply(kins, .ud_probe, grid = grid)
  if (all(vapply(probes, `[[`, logical(1), "safe"))) return(grid)

  env <- .ud_envelope_union(probes)
  box <- .extcalc(pts, ext = ext_factor)
  crs_wkt <- sf::st_crs(pts)$wkt
  cell <- if (!is.null(res)) {
    res
  } else {
    max(box[["xmax"]] - box[["xmin"]], box[["ymax"]] - box[["ymin"]]) / dim_size
  }
  if (is.null(res)) {
    ## fixed point of cell = longest side of the needed range / dim_size.
    ## The needed range grows by at most 2 cells per axis (0.5 + 1.5 in y),
    ## so each step shrinks the gap by a factor of at most 2 / dim_size: it
    ## converges for dim_size >= 3.
    converged <- FALSE
    for (i in seq_len(500)) {
      need <- .ud_needed_range(box, env, cell)
      cell_new <- max(need[["xmax"]] - need[["xmin"]],
                      need[["ymax"]] - need[["ymin"]]) / dim_size
      if (cell_new <= cell * (1 + 1e-12)) {
        converged <- TRUE
        break
      }
      cell <- cell_new
    }
    if (!converged) {
      rlang::abort(sprintf(paste0(
        "`dim_size = %d` is too small to hold the bridges' kernel windows ",
        "with the margin the kernel needs (at least 3). Use a larger ",
        "`dim_size`, or give the cell size as `res`."), as.integer(dim_size)),
        class = "move2utils_ud_dim_size_too_small")
    }
  }
  for (attempt in 0:3) {
    ## A margin of exactly the minimum can lose to floating-point rounding.
    ## With `res` the cell size is fixed, so widen by whole cells; with
    ## `dim_size` the count is fixed, so enlarge the cell a hair instead.
    if (is.null(res) && attempt > 0) cell <- cell * (1 + 1e-9)
    need <- .ud_needed_range(box, env, cell,
                             extra = if (is.null(res)) 0 else attempt)
    if (!is.null(res)) {
      ## With `dim_size` the cell count is fixed; with `res` an enlarged grid
      ## can grow without bound, so cap it before anything is allocated.
      nc <- max(1, ceiling((need[["xmax"]] - need[["xmin"]]) / cell - 1e-9))
      nr <- max(1, ceiling((need[["ymax"]] - need[["ymin"]]) / cell - 1e-9))
      max_cells <- .ud_max_cells()
      if (nc * nr > max_cells) {
        rlang::abort(sprintf(paste0(
          "Enlarging the grid so that every bridge's kernel window fits would ",
          "need %s x %s = %.3g cells of %g map units, above the limit of ",
          "%.3g cells (option `move2utils.ud_max_cells`). Use a larger `res`, ",
          "or give `dim_size`, which fixes the number of cells; or raise the ",
          "limit with options(move2utils.ud_max_cells = ...)."),
          format(nc, scientific = FALSE), format(nr, scientific = FALSE),
          nc * nr, cell, max_cells),
          class = "move2utils_ud_grid_too_large")
      }
    }
    grid <- .raster_from_range(need, cell, crs_wkt)
    if (all(vapply(kins, function(k) .ud_probe(k, grid)$safe, logical(1)))) {
      if (isTRUE(verbose)) {
        rlang::inform(sprintf(paste0(
          "Grid enlarged so that every bridge's kernel window fits: %d x %d ",
          "cells of %.4g map units (the box from ext_factor = %g was too ",
          "small)."), terra::ncol(grid), terra::nrow(grid),
          terra::res(grid)[1], ext_factor),
          class = "move2utils_ud_grid_enlarged")
      }
      return(grid)
    }
  }
  rlang::abort(paste0(
    "Internal error: could not build a UD grid large enough for the kernel ",
    "windows. Please report this with the data and call that produced it."),
    class = "move2utils_ud_grid_internal")
}

#' The largest number of cells an automatic grid may have
#'
#' `getOption("move2utils.ud_max_cells")`, default 1e8 (about 800 MB per
#' layer of doubles).  Every grid built from `res` or `dim_size` is capped,
#' before and after any enlargement; templates are the caller's and are not.
#' @keywords internal
.ud_max_cells <- function() {
  m <- getOption("move2utils.ud_max_cells", 1e8)
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m < 4) {
    rlang::abort(paste0(
      "Option `move2utils.ud_max_cells` must be a single number of at least ",
      "4 (the kernels need two rows and two columns)."),
      class = "move2utils_ud_bad_option")
  }
  m
}

#' Stop before a grid from `res` or `dim_size` exceeds the cell cap
#' @keywords internal
.ud_check_grid_cells <- function(grid, res, dim_size) {
  max_cells <- .ud_max_cells()
  n <- terra::ncell(grid)
  if (n > max_cells) {
    rlang::abort(sprintf(paste0(
      "The grid from %s would have %s x %s = %.3g cells, above the limit ",
      "of %.3g cells (option `move2utils.ud_max_cells`). %s, or raise the ",
      "limit with options(move2utils.ud_max_cells = ...)."),
      if (!is.null(res)) sprintf("`res = %g`", res) else
        sprintf("`dim_size = %d`", as.integer(dim_size)),
      format(terra::ncol(grid), scientific = FALSE),
      format(terra::nrow(grid), scientific = FALSE), n, max_cells,
      if (!is.null(res)) "Use a larger `res`" else "Use a smaller `dim_size`"),
      class = "move2utils_ud_grid_too_large")
  }
  invisible(grid)
}

#' Stop with a useful message when a template raster is too small
#' @keywords internal
.ud_check_template <- function(grid, kins, reprojected) {
  probes <- lapply(kins, .ud_probe, grid = grid)
  if (all(vapply(probes, `[[`, logical(1), "safe"))) return(invisible(TRUE))
  sides <- Reduce(`|`, lapply(probes, `[[`, "sides"))
  side_txt <- paste(names(sides)[sides], collapse = ", ")
  if (terra::nrow(grid) < 2L || terra::ncol(grid) < 2L) {
    rlang::abort(paste0(
      "The raster supplied as `raster`",
      if (reprojected) ", reprojected into the metric CRS used for the computation," else "",
      " has ", terra::nrow(grid), " row(s) and ", terra::ncol(grid),
      " column(s); the UD kernels need at least two of each."),
      class = "move2utils_ud_raster_too_small")
  }
  if (reprojected) {
    rlang::abort(paste0(
      "The raster supplied as `raster`, reprojected into the metric CRS used ",
      "for the computation, is too small: the bridges' kernel windows reach ",
      "beyond its ", side_txt, " edge(s). A template is never enlarged. ",
      "Enlarge it on those sides, or omit `raster` and give `res` or ",
      "`dim_size` to let a large enough grid be built."),
      class = "move2utils_ud_raster_too_small")
  }
  e <- as.vector(terra::ext(grid))
  r <- terra::res(grid)
  need <- .ud_needed_range(c(xmin = e[[1]], xmax = e[[2]], ymin = e[[3]],
                             ymax = e[[4]]),
                           .ud_envelope_union(probes), cell = max(r))
  rlang::abort(sprintf(paste0(
    "The raster supplied as `raster` is too small: the bridges' kernel ",
    "windows reach beyond its %s edge(s). A template is never enlarged. ",
    "This keeps its resolution and alignment and is large enough:\n",
    "  raster <- terra::extend(raster, terra::ext(%s, %s, %s, %s), ",
    "snap = \"out\")\n",
    "Or omit `raster` and give `res` or `dim_size` to let a large enough grid ",
    "be built."), side_txt,
    format(floor(need[["xmin"]]), scientific = FALSE),
    format(ceiling(need[["xmax"]]), scientific = FALSE),
    format(floor(need[["ymin"]]), scientific = FALSE),
    format(ceiling(need[["ymax"]]), scientific = FALSE)),
    class = "move2utils_ud_raster_too_small")
}

#' Canonicalise a move2 object to its per-track local AEQD (metric) CRS,
#' regardless of the input CRS.  This is the stratified projection rule for
#' the outlier functions: geometry is always computed in the same local
#' metric space, so results are invariant to the CRS the caller supplied
#' (longitude/latitude, UTM, another AEQD).  `mt_aeqd_crs()` derives the same
#' centre from the same data in any input CRS, so the transform is a no-op
#' when the object is already in that canonical AEQD (e.g. when a detector is
#' called from a cascade that already projected).  Callers attach flag
#' columns back onto the untouched original object to return the caller's CRS.
#' @keywords internal
.to_canonical_aeqd <- function(x) {
  if (!mt_is_move2(x)) return(x)
  target <- move2::mt_aeqd_crs(x, center = "center", units = "m")
  if (sf::st_crs(x) == target) return(x)
  sf::st_transform(x, target)
}

#' Split a multi-track move2 object into a named list of single-track objects
#' @keywords internal
.split_tracks <- function(x) {
  .validate_move2(x)
  ids <- unique(mt_track_id(x))
  tracks <- lapply(ids, function(id) {
    x[mt_track_id(x) == id, ]
  })
  names(tracks) <- as.character(ids)

  # Validate each track individually
  for (nm in names(tracks)) {
    trk <- tracks[[nm]]
    if (nrow(trk) == 0) next
    coords <- st_coordinates(trk)
    if (any(!is.finite(coords))) {
      rlang::abort(
        sprintf("Track '%s' contains non-finite coordinates.", nm),
        class = "move2utils_input_nonfinite_coords")
    }
    ts_numeric <- as.numeric(mt_time(trk))
    if (is.unsorted(ts_numeric, strictly = FALSE)) {
      rlang::abort(
        sprintf("Track '%s' has timestamps not in chronological order.", nm),
        class = "move2utils_input_unsorted_timestamps")
    }
    if (any(diff(ts_numeric) == 0)) {
      rlang::abort(
        sprintf("Track '%s' contains duplicate timestamps. Remove duplicates first: move2::mt_filter_unique(x)",
                nm),
        class = "move2utils_input_duplicate_timestamps")
    }
  }

  tracks
}

#' Calculate extended bounding box for raster creation
#' @keywords internal
.extcalc <- function(x, ext = 0.3) {
  bb <- st_bbox(x)
  x_range <- as.numeric(bb["xmax"] - bb["xmin"])
  y_range <- as.numeric(bb["ymax"] - bb["ymin"])
  c(
    xmin = as.numeric(bb["xmin"]) - ext * x_range,
    xmax = as.numeric(bb["xmax"]) + ext * x_range,
    ymin = as.numeric(bb["ymin"]) - ext * y_range,
    ymax = as.numeric(bb["ymax"]) + ext * y_range
  )
}

#' Create a raster grid for UD computation
#'
#' @param x A move2 or sf object to derive the extent from.
#' @param cell_size Numeric cell size in map units, or NULL to auto-compute.
#' @param dim_size Number of cells along the longest dimension (used if cell_size is NULL).
#' @param ext Extension factor for the bounding box.
#' @return A `terra::SpatRaster` with the appropriate extent, resolution, and CRS.
#' @keywords internal
.make_raster <- function(x, cell_size = NULL, dim_size = 10, ext = 0.3) {
  range <- .extcalc(x, ext = ext)
  x_range <- range["xmax"] - range["xmin"]
  y_range <- range["ymax"] - range["ymin"]

  if (is.null(cell_size)) {
    cell_size <- max(x_range, y_range) / dim_size
  }

  ymin <- range["ymin"] - (ceiling(y_range / cell_size) * cell_size - y_range) / 2
  ymax <- range["ymax"] + (ceiling(y_range / cell_size) * cell_size - y_range) / 2
  xmin <- range["xmin"] - (ceiling(x_range / cell_size) * cell_size - x_range) / 2
  xmax <- range["xmax"] + (ceiling(x_range / cell_size) * cell_size - x_range) / 2

  nr <- round((ymax - ymin) / cell_size)
  nc <- round((xmax - xmin) / cell_size)

  r <- rast(
    nrows = nr, ncols = nc,
    xmin = xmin, xmax = xmax,
    ymin = ymin, ymax = ymax,
    crs = st_crs(x)$wkt
  )
  r
}


## Fast step-length computation matching mt_distance() convention.
##
## Returns a length-n numeric vector: element i = distance (m) from
## fix i to fix i+1; element n = NA.  On lon/lat input, uses a
## vectorised Haversine with WGS84 mean radius (6371008.8 m); the
## result agrees with `move2::mt_distance()` within ~2 m on continent-
## scale steps and ~mm on typical movement steps, while running about
## 50x faster on 10^5--10^6 location tracks.  On projected input,
## uses direct Euclidean distance from `sf::st_coordinates()`.
##
## NA propagation matches mt_distance(): a non-finite coordinate at
## either end of a step yields NA, and consecutive fixes belonging to
## different tracks yield NA at the cross-track step.
##
## @param x A move2 object (single- or multi-track).
## @return Numeric vector of length nrow(x), with the last element NA.
## @keywords internal
.step_lengths_fast <- function(x) {
  cc <- sf::st_coordinates(x)
  n  <- nrow(cc)
  if (n < 2L) return(rep(NA_real_, n))

  step <- .step_lengths_from_cc(cc, isTRUE(sf::st_is_longlat(x)))

  ## NA at cross-track boundaries when input is multi-track (the raw-
  ## matrix helper has no track-id awareness; we handle that here).
  ids <- move2::mt_track_id(x)
  if (length(unique(ids)) > 1L) {
    ids_chr <- as.character(ids)
    cross   <- ids_chr[-n] != ids_chr[-1L]
    step[c(cross, FALSE)] <- NA_real_
  }
  step
}


## Raw-matrix step-length helper used by `.step_lengths_fast()` and by
## the primitive `.fn_core` entry points that consume `(cc, ...)`
## directly without a `move2` object in scope.
##
## Computes consecutive-pair distances on a 2-column coordinate matrix.
## On `was_longlat = TRUE` uses the WGS84-mean-radius Haversine; on
## `FALSE` uses Cartesian Euclidean.  Returns length-`nrow(cc)`; the
## last element is NA per the package's step-length convention.
## Non-finite coords propagate to NA at both adjacent steps.  Caller is
## responsible for cross-track-boundary NA (the helper assumes single-
## track input).
##
## @keywords internal
.step_lengths_from_cc <- function(cc, was_longlat) {
  n <- nrow(cc)
  if (n < 2L) return(rep(NA_real_, n))

  if (was_longlat) {
    R    <- 6371008.8
    rad  <- pi / 180
    lon1 <- cc[-n,   1] * rad
    lat1 <- cc[-n,   2] * rad
    lon2 <- cc[-1L,  1] * rad
    lat2 <- cc[-1L,  2] * rad
    dlat <- lat2 - lat1
    dlon <- lon2 - lon1
    a    <- sin(dlat / 2) ^ 2 +
            cos(lat1) * cos(lat2) * sin(dlon / 2) ^ 2
    step <- 2 * R * asin(pmin(sqrt(a), 1))
  } else {
    step <- sqrt(diff(cc[, 1]) ^ 2 + diff(cc[, 2]) ^ 2)
  }

  bad <- !is.finite(cc[-n, 1]) | !is.finite(cc[-n, 2]) |
         !is.finite(cc[-1L, 1]) | !is.finite(cc[-1L, 2])
  step[bad] <- NA_real_
  c(step, NA_real_)
}


## Vectorised great-circle distance (m) between paired lon/lat
## points.  Same WGS84 mean radius (6371008.8 m) and clamp pattern
## as `.step_lengths_fast`; intended for non-adjacent pair distances
## such as the (i-k, i+k) chord in `mt_flag_outliers_detour()`.
##
## Inputs are decimal-degree numeric vectors of equal length.
## Output is the same length.  Non-finite inputs propagate to NA.
##
## Choice of Haversine over project-then-Euclidean for non-adjacent
## chords on lon/lat: at 10^6+ fixes the AEQD projection cost is
## the dominant single line of work, while a vectorised Haversine
## pair distance runs in tens of nanoseconds per pair.  Haversine
## also avoids AEQD distortion on continental-scale tracks where
## the centroid-centred projection becomes inaccurate at the
## periphery.
##
## @keywords internal
.haversine_pair <- function(lon1, lat1, lon2, lat2) {
  R    <- 6371008.8
  rad  <- pi / 180
  phi1 <- lat1 * rad
  phi2 <- lat2 * rad
  dphi <- (lat2 - lat1) * rad
  dlmb <- (lon2 - lon1) * rad
  a    <- sin(dphi / 2) ^ 2 +
          cos(phi1) * cos(phi2) * sin(dlmb / 2) ^ 2
  d    <- 2 * R * asin(pmin(sqrt(a), 1))
  bad  <- !is.finite(lon1) | !is.finite(lat1) |
          !is.finite(lon2) | !is.finite(lat2)
  d[bad] <- NA_real_
  d
}


## Fast turn-angle computation matching move2::mt_turnangle() convention.
##
## Returns a length-n numeric vector in radians: element 1 = NA;
## element i (2 <= i <= n-1) = turn angle at fix i (the change in
## heading between the incoming step (i-1 -> i) and outgoing step
## (i -> i+1)); element n = NA.  Wrapped to [-pi, pi].
##
## On lon/lat input, computes great-circle azimuth using the standard
## spherical formula.  This differs from move2::mt_turnangle(), which
## uses the WGS84 ellipsoid via s2, by < 0.003 rad in the worst case
## (typically < 0.001 rad), well below the precision that any
## downstream KDE-based scoring is sensitive to.  In return, this
## helper runs about 25x faster on 10^5--10^6 location tracks.  On
## projected input, computes Cartesian atan2(dx, dy) -- works directly
## without the move2 limitation that mt_azimuth() rejects projected
## coordinates.
##
## NA propagation: any cross-track step yields NA at both adjacent
## fixes; non-finite coordinates propagate similarly.
##
## @param x A move2 object (single- or multi-track).
## @return Numeric vector of length nrow(x), turn angles in radians.
## @keywords internal
.turn_angles_fast <- function(x) {
  cc <- sf::st_coordinates(x)
  turn <- .turn_angles_from_cc(cc, isTRUE(sf::st_is_longlat(x)))

  ## NA at cross-track boundaries when input is multi-track
  ids <- move2::mt_track_id(x)
  if (length(unique(ids)) > 1L) {
    n <- nrow(cc)
    ids_chr <- as.character(ids)
    cross   <- ids_chr[-n] != ids_chr[-1L]
    bad_in  <- c(NA, cross)
    bad_out <- c(cross, NA)
    turn[bad_in | bad_out] <- NA_real_
  }
  turn
}


## Raw-matrix turn-angle helper used by `.turn_angles_fast()` and by
## primitive `.fn_core` entry points that consume `(cc, ...)` directly.
## Same convention and accuracy as `.turn_angles_fast()`; assumes
## single-track input (caller handles cross-track NA propagation).
##
## @keywords internal
.turn_angles_from_cc <- function(cc, was_longlat) {
  n  <- nrow(cc)
  if (n < 3L) return(rep(NA_real_, n))

  if (was_longlat) {
    rad  <- pi / 180
    lon  <- cc[, 1] * rad
    lat  <- cc[, 2] * rad
    lon1 <- lon[-n];  lat1 <- lat[-n]
    lon2 <- lon[-1L]; lat2 <- lat[-1L]
    dlon <- lon2 - lon1
    y_a  <- sin(dlon) * cos(lat2)
    x_a  <- cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon)
    az   <- atan2(y_a, x_a)
  } else {
    az <- atan2(diff(cc[, 1]), diff(cc[, 2]))
  }

  turn <- rep(NA_real_, n)
  turn[2:(n - 1L)] <- az[2:(n - 1L)] - az[1:(n - 2L)]
  turn <- ((turn + pi) %% (2 * pi)) - pi

  bad_coord <- !is.finite(cc[, 1]) | !is.finite(cc[, 2])
  if (any(bad_coord)) {
    bad_in  <- c(NA, bad_coord[-n] | bad_coord[-1L])
    bad_out <- c(bad_coord[-n] | bad_coord[-1L], NA)
    turn[bad_in | bad_out] <- NA_real_
  }

  turn
}

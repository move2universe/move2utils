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
#' @keywords internal
.ud_resolve_grid <- function(raster, td, target_crs, restore, dim_size, ext) {
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
      return(list(grid = terra::project(raster, compute_crs$wkt),
                  warp_to = raster))
    }
    warp_to <- if (isTRUE(restore) && target != compute_crs) target else NULL
    return(list(grid = raster, warp_to = warp_to))
  }

  ## auto grid (NULL or numeric cell size) built in the metric compute CRS
  pts <- sf::st_as_sf(data.frame(x = td$x, y = td$y),
                      coords = c("x", "y"), crs = compute_crs)
  grid <- if (is.numeric(raster) && length(raster) == 1L) {
    .make_raster(pts, cell_size = raster, ext = ext)
  } else {
    .make_raster(pts, dim_size = dim_size, ext = ext)
  }
  warp_to <- if (isTRUE(restore) && target != compute_crs) target else NULL
  list(grid = grid, warp_to = warp_to)
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

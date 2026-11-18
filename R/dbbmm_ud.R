#' Dynamic Brownian Bridge Movement Model — Utilisation Distribution
#'
#' Compute a utilisation distribution (UD) from a movement track using
#' the dynamic Brownian bridge movement model. Accepts a `move2` object
#' directly (estimates variance internally), a pre-computed variance
#' object, or a named list of variance objects for multi-track input.
#'
#' @param object One of:
#'   \itemize{
#'     \item A `move2` object (single or multi-track; any CRS — see
#'       \dQuote{Coordinate reference systems} below). Variance is
#'       estimated internally via [mt_dbbmm_variance()].
#'     \item An `mt_dbbmm_variance` object from [mt_dbbmm_variance()].
#'     \item A named list of `mt_dbbmm_variance` objects (multi-track).
#'   }
#' @param raster A `terra::SpatRaster` template: the UD is computed for its
#'   cells and returned on that grid. A template is never enlarged; if it is
#'   too small for the bridges, the call stops and says on which sides; when
#'   the template is in the CRS the computation runs in, it also gives the
#'   `terra::extend()` call that fixes it. Give at most one of `raster`,
#'   `res` and `dim_size`; with none, a grid is built from `dim_size = 100`.
#'   (A number here is the deprecated way of giving `res`.)
#' @param res Numeric. Cell size of an automatically built grid, in map
#'   units of the compute CRS (metres). Cells are square. For multi-track
#'   input, one common grid covers all tracks. A grid built from `res` or
#'   `dim_size`, enlarged or not, may have at most
#'   `getOption("move2utils.ud_max_cells", 1e8)` cells; above that the call
#'   stops before anything is allocated (class
#'   `move2utils_ud_grid_too_large`) and suggests a larger `res`.
#' @param location_error Per-fix horizontal 1-sigma positional error,
#'   **in metres**.  When \code{object} is a \code{move2} object, this is
#'   forwarded to [mt_dbbmm_variance()] (see that function for the full
#'   set of accepted forms: NULL / scalar / vector / column name /
#'   \code{"auto"}).  When \code{object} is an \code{mt_dbbmm_variance}
#'   object or a list of them, the default \code{NULL} re-uses the
#'   per-fix vector stored on each variance object at fit time;
#'   supplying an explicit value here overrides it.
#' @param location_error_na Strategy for filling \code{NA}s in a per-fix
#'   \code{location_error} vector.  Forwarded to [mt_dbbmm_variance()]
#'   when \code{object} is a \code{move2}; ignored otherwise (the
#'   variance object already carries the imputed vector).
#' @param margin Integer (odd). Margin for variance estimation window.
#'   Only used when `object` is a `move2` object.
#' @param window_size Integer (odd). Window size for variance estimation.
#'   Only used when `object` is a `move2` object.
#' @param ext_factor Numeric, non-negative. Fraction of the track's
#'   bounding-box width and height added on each side of an automatically
#'   built grid. Default 0.5. If that box is too small for the bridges, the
#'   grid is enlarged to fit them (see \dQuote{Grid}).
#' @param dim_size Integer, at least 2. Number of cells along the longest
#'   side of an automatically built grid. Higher values give finer
#'   resolution but slower computation. Used when neither `raster` nor
#'   `res` is given, with a default of 100.
#' @param ext Deprecated former name of `ext_factor`. Renamed because in
#'   `terra` `ext` denotes an extent.
#' @param time_step Numeric or `NULL`. Time step for the Brownian bridge
#'   integration, in minutes. Defaults to 1/15 of the minimum time lag.
#'   Smaller values give more precise UDs but are slower.
#' @param verbose Logical. If `TRUE`, print a computational size estimate,
#'   report progress for multi-track input, and say when an automatic grid
#'   had to be enlarged.
#' @param ... Additional arguments passed to methods.
#'
#' @return For single-track input: a `terra::SpatRaster` where cell
#'   values sum to 1.0, representing the utilisation distribution.
#'
#'   For multi-track input: a multi-layer `terra::SpatRaster` with one
#'   named layer per track, all on a common grid. Each layer sums to 1.0.
#'
#' @details
#' The UD is computed by evaluating the Brownian bridge probability
#' density at regular time steps along each segment and accumulating
#' the density onto a raster grid. The computation is implemented in C
#' with OpenMP parallelisation of the inner grid loops for performance.
#'
#' For multi-track input, a common raster grid is computed from the
#' combined spatial extent of all tracks, ensuring that UDs are
#' directly comparable.
#'
#' @section Coordinate reference systems:
#' The UD is returned in the same CRS you provide, so it overlays directly
#' on your other layers: give longitude/latitude and you get a long/lat UD,
#' give UTM and you get a UTM UD. The Brownian-bridge math needs metric
#' coordinates, so when the target CRS is already projected in metres it is
#' used directly (lossless); longitude/latitude (or any non-metric CRS) is
#' computed in a local AEQD and the resulting UD is then reprojected back to
#' your CRS and renormalised so it still sums to 1. If you pass an
#' environmental raster as `raster`, the UD is returned on \emph{that} grid
#' (reprojecting the movement internally if needed), ready to stack against
#' your environmental layers. The reprojection of a long/lat UD resamples
#' the surface; supply your data in a suitable metric projection if you want
#' to avoid that step entirely.
#'
#' @section Grid:
#' The grid is either the template given as `raster`, or it is built
#' automatically: the track's bounding box, widened by `ext_factor`, divided
#' into square cells of size `res` or into `dim_size` cells along its longest
#' side. The kernel spreads each bridge time step over four standard
#' deviations around the bridge mean, so the grid has to reach that far
#' beyond the track. An automatic grid that does not is enlarged until it
#' does: with `res` the cell size holds and cells are added; with `dim_size`
#' the cell count holds and the cells grow. A grid that is large enough is
#' used unchanged, so results are identical to earlier versions wherever
#' those succeeded. A template is never enlarged: if it is too small the
#' call stops and says on which sides. Grids need at least two rows and two
#' columns.
#'
#' @section Coming from move:
#' `mt_dbbmm_ud()` is the port of `move::brownian.bridge.dyn()`, and its C
#' kernel is `move`'s. Arguments map as follows: a `RasterLayer` as `raster`
#' becomes a `SpatRaster` template as `raster`; a number as `raster` (the
#' cell size) becomes `res`; `dimSize` becomes `dim_size`; `ext` becomes
#' `ext_factor`; `bbox` becomes a template,
#' `raster = terra::rast(<extent>, resolution = <cell size>, crs = <crs>)`.
#' The defaults differ: `move` uses `ext = 0.3` and `dimSize = 10`, this
#' package `ext_factor = 0.5` and `dim_size = 100`. Where `move` stops
#' because the grid is too small, an automatic grid is enlarged here
#' instead. `move`'s kernel also accepts a grid on which its window reaches
#' exactly one cell beyond the last row or column, and then reads and writes
#' outside the grid; such a grid counts as too small here. Given the same
#' grid and the same variance, the UD agrees with `move` to within 3e-9
#' relative (measured on `move::leroy`).
#'
#' @references
#' Kranstauber, B., Kays, R., LaPoint, S. D., Wikelski, M., & Safi, K.
#' (2012). A dynamic Brownian bridge movement model to estimate utilization
#' distributions for heterogeneous animal movement. *Journal of Animal
#' Ecology*, 81(4), 738-746. \doi{10.1111/j.1365-2656.2012.01955.x}
#'
#' @seealso [mt_dbbmm_variance()] to compute variance separately,
#'   [mt_dbgb_ud()] for the directional variant, [mt_motion_variance()]
#'   to extract variances.
#'
#' @examples
#' library(move2)
#' library(move2utils)
#' library(sf)
#' fishers <- mt_read(mt_example())
#' fishers <- fishers[!st_is_empty(fishers), ]
#' leroy <- filter_track_data(fishers, .track_id = "M1")
#' ## the first 200 fixes keep the example quick
#' leroy <- st_transform(leroy, mt_aeqd_crs(leroy))[seq_len(200), ]
#' ## One step: variance and UD
#' ud <- mt_dbbmm_ud(leroy, location_error = 25,
#'                   window_size = 31, margin = 11)
#' terra::plot(ud)
#'
#' ## Two steps: fit the variance once, re-use it for other grids
#' var_obj <- mt_dbbmm_variance(leroy, location_error = 25,
#'                              window_size = 31, margin = 11)
#' ud_fine <- mt_dbbmm_ud(var_obj, dim_size = 300, ext_factor = 1.25)
#' ud_50m  <- mt_dbbmm_ud(var_obj, res = 50)
#'
#' @export
mt_dbbmm_ud <- function(object,
                         raster = NULL,
                         location_error = NULL,
                         margin = 11,
                         window_size = 31,
                         ext_factor = 0.5,
                         dim_size = NULL,
                         time_step = NULL,
                         verbose = TRUE,
                         location_error_na = "median",
                         res = NULL,
                         ext = NULL,
                         ...) {
  UseMethod("mt_dbbmm_ud")
}

#' @export
mt_dbbmm_ud.move2 <- function(object,
                                raster = NULL,
                                location_error = NULL,
                                margin = 11,
                                window_size = 31,
                                ext_factor = 0.5,
                                dim_size = NULL,
                                time_step = NULL,
                                verbose = TRUE,
                                location_error_na = "median",
                                res = NULL,
                                ext = NULL,
                                ...) {
  ga <- .ud_grid_args(raster, res, dim_size, ext_factor, ext,
                      ext_factor_given = !missing(ext_factor))
  raster <- ga$raster
  ## CRS stratification (see ?mt_dbbmm_ud "Coordinate reference systems"):
  ## the UD is returned in the TARGET crs -- the env template's crs if one was
  ## supplied, else the movement's own input crs.  The bridge math runs in a
  ## metric COMPUTE crs: the target itself when it is suitable (projected,
  ## metres), else a local AEQD.  Projecting movement onto the compute crs
  ## here (rather than inside the variance call) keeps the suitable-target
  ## fast path lossless: compute == target, so no warp at the end.
  target <- if (inherits(raster, "SpatRaster")) {
    sf::st_crs(terra::crs(raster))
  } else {
    sf::st_crs(object)
  }
  compute_crs <- if (.crs_is_suitable(target)) {
    target
  } else {
    sf::st_crs(move2::mt_aeqd_crs(object, center = "center", units = "m"))
  }
  if (sf::st_crs(object) != compute_crs) {
    object <- sf::st_transform(object, compute_crs)
  }
  var_obj <- mt_dbbmm_variance(object, location_error = location_error,
                                window_size = window_size, margin = margin,
                                location_error_na = location_error_na)
  ## Pass location_error = NULL so the variance object's stored per-fix
  ## vector is used unchanged.  Avoids re-passing a full-length user
  ## vector that no longer matches the per-track variance n_locs.
  mt_dbbmm_ud(var_obj, raster = ga$raster, res = ga$res,
               dim_size = ga$dim_size, location_error = NULL,
               ext_factor = ga$ext_factor, time_step = time_step,
               verbose = verbose, .target_crs = target)
}

#' @export
mt_dbbmm_ud.list <- function(object,
                               raster = NULL,
                               location_error = NULL,
                               margin = 11,
                               window_size = 31,
                               ext_factor = 0.5,
                               dim_size = NULL,
                               time_step = NULL,
                               verbose = TRUE,
                               location_error_na = "median",
                               res = NULL,
                               ext = NULL,
                               ...,
                               .target_crs = NULL) {
  # Validate: must be a named list of mt_dbbmm_variance objects
  if (!all(vapply(object, inherits, logical(1), "mt_dbbmm_variance"))) {
    rlang::abort("List must contain mt_dbbmm_variance objects.",
                 class = "move2utils_dbbmm_ud_bad_list")
  }
  ga <- .ud_grid_args(raster, res, dim_size, ext_factor, ext,
                      ext_factor_given = !missing(ext_factor))

  ## location_error handling:
  ##   NULL                  -> each variance object uses its stored vector
  ##   numeric scalar / NULL -> applied uniformly
  ##   named list per track  -> per-track override
  ## Per-row vectors are intentionally rejected here -- the variance
  ## objects carry no row-mapping back to a flat user vector, so the
  ## slice would be ambiguous.  Build the variance list with per-track
  ## location_error if you need per-fix overrides.
  if (!is.null(location_error)) {
    if (is.numeric(location_error) && length(location_error) == 1L) {
      le_per_track <- setNames(as.list(rep(location_error, length(object))),
                                names(object))
    } else if (is.list(location_error) &&
               !is.null(names(location_error)) &&
               all(names(object) %in% names(location_error))) {
      le_per_track <- location_error[names(object)]
    } else {
      rlang::abort(paste0(
        "With a list of variance objects, `location_error` must be NULL ",
        "(use stored per-track vectors), a numeric scalar (uniform across ",
        "tracks), or a named list keyed by track id. Per-row vectors ",
        "aren't accepted at the list-dispatch step -- pass them to ",
        "mt_dbbmm_variance() instead."),
        class = "move2utils_dbbmm_ud_bad_location_error")
    }
  } else {
    le_per_track <- setNames(vector("list", length(object)), names(object))
  }

  ## CRS stratification: all tracks share one metric COMPUTE crs (the
  ## variance objects were built from a single multi-track object, so
  ## .ensure_projected centred one AEQD on the whole object).  Build the
  ## common grid and compute every layer in that compute crs, then warp the
  ## whole stack ONCE to the target crs -- warping per layer would let the
  ## reprojected grids diverge and break the stack.
  compute_crs <- sf::st_crs(object[[1]]$track_data$crs)
  user_template <- ga$raster

  ## the common grid has to hold every track's kernel windows
  kins <- lapply(names(object), function(nm) {
    .dbbmm_kernel_inputs(object[[nm]], le_per_track[[nm]], time_step)
  })

  if (is.null(user_template)) {
    all_x <- unlist(lapply(object, function(v) v$track_data$x))
    all_y <- unlist(lapply(object, function(v) v$track_data$y))
    pts <- sf::st_as_sf(data.frame(x = all_x, y = all_y),
                         coords = c("x", "y"), crs = compute_crs)
    common_raster <- .ud_auto_grid(pts, ga, kins, verbose = verbose)
  } else if (sf::st_crs(terra::crs(user_template)) != compute_crs) {
    ## explicit env grid in another crs -> reproject into the compute crs for
    ## the kernel; the stack is warped back onto the env grid below
    common_raster <- terra::project(user_template, compute_crs$wkt)
    .ud_check_template(common_raster, kins, reprojected = TRUE)
  } else {
    common_raster <- user_template
    .ud_check_template(common_raster, kins, reprojected = FALSE)
  }

  # Compute UD for each track on the common grid, in the compute crs
  layers <- lapply(names(object), function(nm) {
    if (verbose) message("Computing UD for track: ", nm)
    mt_dbbmm_ud(object[[nm]], raster = common_raster,
                 location_error = le_per_track[[nm]],
                 ext_factor = ga$ext_factor,
                 time_step = time_step, verbose = FALSE,
                 .restore = FALSE)
  })

  # Stack into multi-layer SpatRaster, then warp the stack to the target crs
  stk <- terra::rast(layers)
  names(stk) <- names(object)

  target <- if (!is.null(.target_crs)) {
    sf::st_crs(.target_crs)
  } else if (!is.null(user_template)) {
    sf::st_crs(terra::crs(user_template))
  } else if (!is.null(object[[1]]$track_data$orig_crs)) {
    sf::st_crs(object[[1]]$track_data$orig_crs)
  } else {
    compute_crs
  }
  if (target != compute_crs) {
    warp_to <- if (!is.null(user_template)) user_template else target
    stk <- .warp_ud_to_target(stk, warp_to)
    names(stk) <- names(object)
  }
  stk
}

#' @export
mt_dbbmm_ud.mt_dbbmm_variance <- function(object,
                                            raster = NULL,
                                            location_error = NULL,
                                            margin = 11,
                                            window_size = 31,
                                            ext_factor = 0.5,
                                            dim_size = NULL,
                                            time_step = NULL,
                                            verbose = TRUE,
                                            location_error_na = "median",
                                            res = NULL,
                                            ext = NULL,
                                            ...,
                                            .target_crs = NULL,
                                            .restore = TRUE) {
  ga <- .ud_grid_args(raster, res, dim_size, ext_factor, ext,
                      ext_factor_given = !missing(ext_factor))
  td <- object$track_data

  ## The kernel's exact inputs, shared with the grid-bounds probe so that the
  ## grid is checked against precisely what the kernel will evaluate.
  kin <- .dbbmm_kernel_inputs(object, location_error, time_step)

  ## CRS stratification: build the kernel grid in the metric compute crs
  ## (td$crs) and learn where to warp the result (see .ud_resolve_grid).
  grid_res <- .ud_resolve_grid(ga, list(kin), td, .target_crs, .restore,
                               verbose = verbose)
  raster   <- grid_res$grid

  if (verbose) {
    comp_size <- terra::ncell(raster) *
      (sum(kin$time_lag[object$interest]) / kin$time_step)
    message(sprintf("Computational size: %.1e", comp_size))
  }

  gv <- .ud_grid_vectors(raster, "dbbmm")
  ans <- .Call(
    "dbbmm2_omp",
    kin$x, kin$y, kin$s1, kin$t, kin$le, gv$x, gv$y,
    kin$time_step, kin$ext, kin$interest
  )

  total <- sum(ans)
  if (total == 0 || !is.finite(total)) {
    rlang::abort(paste0(
      "UD computation produced zero or non-finite values. ",
      "The raster extent may not overlap the track."),
      class = "move2utils_dbbmm_ud_zero_total_mass")
  }
  ans <- ans / total
  terra::values(raster) <- ans

  ## return the UD in the target crs (warp + renormalise) when it differs
  ## from the compute crs; otherwise it is already on the caller's grid
  if (!is.null(grid_res$warp_to)) {
    raster <- .warp_ud_to_target(raster, grid_res$warp_to)
  }
  raster
}

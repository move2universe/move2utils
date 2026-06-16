#' Dynamic Bivariate Gaussian Bridge — Utilisation Distribution
#'
#' Compute a utilisation distribution using the dynamic bivariate Gaussian
#' bridge model, which decomposes movement variance into parallel and
#' orthogonal components relative to the direction of travel.
#'
#' @param object One of:
#'   \itemize{
#'     \item A `move2` object (single or multi-track; any CRS — see
#'       \dQuote{Coordinate reference systems} below).
#'     \item An `mt_dbgb_variance` object from [mt_dbgb_variance()].
#'     \item A named list of `mt_dbgb_variance` objects (multi-track).
#'   }
#' @param raster A `terra::SpatRaster` defining the output grid, a numeric
#'   cell size in map units, or `NULL` to auto-compute.
#' @param location_error Per-fix horizontal 1-sigma positional error,
#'   **in metres**.  When \code{object} is a \code{move2} object, this
#'   is forwarded to [mt_dbgb_variance()] (see that function for the
#'   full set of accepted forms: NULL / scalar / vector / column name /
#'   \code{"auto"}).  When \code{object} is an \code{mt_dbgb_variance}
#'   object or a list of them, the default \code{NULL} re-uses the
#'   per-fix vector stored on each variance object at fit time;
#'   supplying an explicit value here overrides it.
#' @param location_error_na Strategy for filling \code{NA}s in a per-fix
#'   \code{location_error} vector.  Forwarded to [mt_dbgb_variance()]
#'   when \code{object} is a \code{move2}; ignored otherwise (the
#'   variance object already carries the imputed vector).
#' @param margin Integer (odd). Margin for variance estimation. Only used
#'   when `object` is a `move2` object.
#' @param window_size Integer (odd). Window size for variance estimation.
#'   Only used when `object` is a `move2` object.
#' @param ext Numeric. Extension factor for the bounding box. Default 0.5.
#' @param dim_size Integer. Cells along the longest dimension. Default 100.
#' @param time_step Numeric or `NULL`. Integration time step in minutes.
#' @param verbose Logical. Print progress messages for multi-track.
#' @param ... Additional arguments passed to methods.
#'
#' @return For single-track: a `terra::SpatRaster` (values sum to 1.0).
#'   For multi-track: a multi-layer `terra::SpatRaster` on a common grid,
#'   one named layer per track, each summing to 1.0.
#'
#' @details
#' The dBGB UD uses an anisotropic Gaussian kernel at each time step,
#' with the kernel elongated along the direction of travel. This produces
#' narrower UDs along directed segments and wider UDs where movement is
#' more random, better capturing the actual space use of the animal.
#'
#' @section Coordinate reference systems:
#' The UD is returned in the CRS you provide (long/lat in → long/lat out,
#' UTM in → UTM out), so it overlays directly on your other layers. The
#' bridge math runs in metres: a projected metric CRS is used directly,
#' while longitude/latitude is computed in a local AEQD and reprojected back
#' (and renormalised to sum 1). Passing an environmental raster as `raster`
#' returns the UD on that exact grid. See [mt_dbbmm_ud()] for the full
#' description.
#'
#' @references
#' Kranstauber, B., Safi, K., & Bartumeus, F. (2014). Bivariate Gaussian
#' bridges: directional factorization of diffusion in Brownian bridge
#' models. *Movement Ecology*, 2(1), 5. \doi{10.1186/2051-3933-2-5}
#'
#' @seealso [mt_dbgb_variance()] to compute variance separately,
#'   [mt_dbbmm_ud()] for the isotropic variant.
#'
#' @examples
#' \dontrun{
#' library(move2)
#' library(sf)
#'
#' fishers <- mt_read(mt_example())
#' fishers <- fishers[!st_is_empty(fishers), ]
#' f1 <- fishers[mt_track_id(fishers) == "F1", ]
#' f1_proj <- st_transform(f1, mt_aeqd_crs(f1))
#'
#' ud <- mt_dbgb_ud(f1_proj, location_error = 25,
#'                   margin = 15, window_size = 31, ext = 0.85)
#' terra::plot(ud)
#' }
#'
#' @export
mt_dbgb_ud <- function(object,
                        raster = NULL,
                        location_error = NULL,
                        margin = 15,
                        window_size = 31,
                        ext = 0.5,
                        dim_size = 100,
                        time_step = NULL,
                        verbose = TRUE,
                        location_error_na = "median",
                        ...) {
  UseMethod("mt_dbgb_ud")
}

#' @export
mt_dbgb_ud.move2 <- function(object,
                               raster = NULL,
                               location_error = NULL,
                               margin = 15,
                               window_size = 31,
                               ext = 0.5,
                               dim_size = 100,
                               time_step = NULL,
                               verbose = TRUE,
                               location_error_na = "median",
                               ...) {
  ## CRS stratification (see ?mt_dbgb_ud "Coordinate reference systems"):
  ## return the UD in the TARGET crs (env template's crs if supplied, else the
  ## movement's input crs); run the bridge math in a metric COMPUTE crs (the
  ## target if suitable, else a local AEQD).  Projecting movement here keeps
  ## the suitable-target path lossless (compute == target -> no warp).
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
  var_obj <- mt_dbgb_variance(object, location_error = location_error,
                               margin = margin, window_size = window_size,
                               location_error_na = location_error_na)
  ## Pass location_error = NULL so each variance object's stored per-fix
  ## vector is used unchanged.
  mt_dbgb_ud(var_obj, raster = raster, location_error = NULL,
              ext = ext, dim_size = dim_size, time_step = time_step,
              verbose = verbose, .target_crs = target)
}

#' @export
mt_dbgb_ud.list <- function(object,
                              raster = NULL,
                              location_error = NULL,
                              margin = 15,
                              window_size = 31,
                              ext = 0.5,
                              dim_size = 100,
                              time_step = NULL,
                              verbose = TRUE,
                              ...,
                              .target_crs = NULL) {
  if (!all(vapply(object, inherits, logical(1), "mt_dbgb_variance"))) {
    rlang::abort("List must contain mt_dbgb_variance objects.",
                 class = "move2utils_dbgb_ud_bad_list")
  }

  ## See mt_dbbmm_ud.list for the location_error contract.  Same rules.
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
        "mt_dbgb_variance() instead."),
        class = "move2utils_dbgb_ud_bad_location_error")
    }
  } else {
    le_per_track <- setNames(vector("list", length(object)), names(object))
  }

  ## CRS stratification: compute every layer in the shared metric COMPUTE crs
  ## on a common grid, then warp the whole stack ONCE to the target (warping
  ## per layer would let the reprojected grids diverge).  See mt_dbbmm_ud.list.
  compute_crs <- sf::st_crs(object[[1]]$track_data$crs)
  user_template <- if (inherits(raster, "SpatRaster")) raster else NULL

  if (is.null(raster) || (is.numeric(raster) && length(raster) == 1)) {
    all_x <- unlist(lapply(object, function(v) v$track_data$x))
    all_y <- unlist(lapply(object, function(v) v$track_data$y))
    pts <- sf::st_as_sf(data.frame(x = all_x, y = all_y),
                         coords = c("x", "y"), crs = compute_crs)
    if (is.numeric(raster)) {
      common_raster <- .make_raster(pts, cell_size = raster, ext = ext)
    } else {
      common_raster <- .make_raster(pts, dim_size = dim_size, ext = ext)
    }
  } else if (sf::st_crs(terra::crs(raster)) != compute_crs) {
    common_raster <- terra::project(raster, compute_crs$wkt)
  } else {
    common_raster <- raster
  }

  layers <- lapply(names(object), function(nm) {
    if (verbose) message("Computing UD for track: ", nm)
    mt_dbgb_ud(object[[nm]], raster = common_raster,
                location_error = le_per_track[[nm]],
                ext = ext, dim_size = dim_size,
                time_step = time_step, verbose = FALSE,
                .restore = FALSE)
  })

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
mt_dbgb_ud.mt_dbgb_variance <- function(object,
                                          raster = NULL,
                                          location_error = NULL,
                                          margin = 15,
                                          window_size = 31,
                                          ext = 0.5,
                                          dim_size = 100,
                                          time_step = NULL,
                                          verbose = TRUE,
                                          ...,
                                          .target_crs = NULL,
                                          .restore = TRUE) {
  td <- object$track_data
  n <- td$n_locs
  if (is.null(location_error)) {
    location_error <- if (!is.null(td$location_error)) {
      td$location_error
    } else {
      rep(0, n)
    }
  } else {
    location_error <- .expand_loc_error(location_error, n)
  }

  # seg_interest is a per-location include mask: seg_interest[i] = TRUE
  # means "include the bridge from location i to location i+1 in the UD".
  # In the natural (auto-computed) case the mask is symmetric around the
  # track midpoint, so previous versions of this routine ORed it with its
  # reverse as a defensive idiom; that was a no-op on natural masks but
  # silently undid user-level masking (see mt_mask_segments()).
  points_interest <- object$seg_interest

  ## CRS stratification: kernel grid in the metric compute crs + warp target.
  grid_res <- .ud_resolve_grid(raster, td, .target_crs, .restore,
                               dim_size, ext)
  raster   <- grid_res$grid

  t_mins <- td$time_mins
  if (is.null(time_step)) {
    time_step <- min(diff(t_mins)) / 20.1
  }

  para_sd <- object$para_sd
  orth_sd <- object$orth_sd
  para_sd[is.na(para_sd)] <- 0
  orth_sd[is.na(orth_sd)] <- 0

  x_grid <- terra::xFromCol(raster, 1:terra::ncol(raster))
  y_grid <- sort(unique(terra::yFromRow(raster, 1:terra::nrow(raster))))

  ans <- .Call(
    "bgb_omp",
    td$x[points_interest], td$y[points_interest],
    para_sd[points_interest], orth_sd[points_interest],
    t_mins[points_interest],
    rep(location_error, length.out = n)[points_interest],
    x_grid, y_grid, time_step, 5
  )

  total <- sum(ans)
  if (total == 0 || !is.finite(total)) {
    rlang::abort(paste0(
      "UD computation produced zero or non-finite values. ",
      "The raster extent may not overlap the track, or ext may be too large."),
      class = "move2utils_dbgb_ud_zero_total_mass")
  }
  ans <- ans / total
  terra::values(raster) <- ans

  if (!is.null(grid_res$warp_to)) {
    raster <- .warp_ud_to_target(raster, grid_res$warp_to)
  }
  raster
}

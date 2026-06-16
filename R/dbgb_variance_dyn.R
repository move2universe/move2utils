#' Dynamic Bivariate Gaussian Bridge Variance Estimation
#'
#' Estimate dynamic parallel and orthogonal movement variances using a
#' sliding window approach with BIC-based breakpoint detection. The
#' bivariate Gaussian bridge decomposes the movement variance into a
#' component parallel to the direction of travel and one orthogonal to
#' it, providing more detail about the movement process than the
#' isotropic dBBMM.
#'
#' @param object A `move2` object. A projected (metric) CRS is used as
#'   supplied; longitude/latitude input is auto-projected to a local AEQD
#'   internally (variance values are CRS-invariant metres either way).
#'   Both single-track and multi-track objects are accepted.
#' @param location_error Per-fix horizontal 1-sigma positional error,
#'   **in metres** (the AEQD / projected-CRS unit). Treated as the
#'   anchor measurement-error prior in the bivariate Brownian-bridge
#'   model.
#'
#'   Accepted forms (identical contract to [mt_dbbmm_variance()]):
#'   \code{NULL}, a non-negative numeric scalar, a numeric vector of
#'   length \code{nrow(object)}, a column name in \code{object}, or
#'   \code{"auto"} (probes \code{eobs_horizontal_accuracy_estimate},
#'   then \code{argos_lc}).  Multi-track inputs are sliced per-track
#'   internally; \code{NA}s are imputed via \code{location_error_na}.
#' @param location_error_na Strategy for filling \code{NA}s in a per-fix
#'   \code{location_error} vector.  One of \code{"median"} (default),
#'   \code{"mean"}, \code{"zero"} (fakes certainty -- use deliberately),
#'   or \code{"approx"} (linear interpolation).  See
#'   [mt_dbbmm_variance()] for the rationale.
#' @param window_size Integer (must be odd). Number of locations in the
#'   sliding window.
#' @param margin Integer (must be odd). Minimum locations on each side of
#'   a potential breakpoint.
#' @param parallel,cores Retained for API compatibility. Since 0.2.0
#'   the variance kernel is implemented in C with OpenMP threading
#'   (commit `bf8eb34`); these arguments are no-ops at the R level.
#'   To control thread count, set the `OMP_NUM_THREADS` environment
#'   variable before invoking R.
#'
#' @return For single-track input: an `mt_dbgb_variance` S3 object
#'   containing:
#'   \describe{
#'     \item{`para_sd`}{Numeric vector of parallel standard deviations.}
#'     \item{`orth_sd`}{Numeric vector of orthogonal standard deviations.}
#'     \item{`n_estim`}{Number of windows each position was estimated in.}
#'     \item{`seg_interest`}{Logical vector of fully-covered segments.}
#'     \item{`margin`, `window_size`}{Parameters used.}
#'     \item{`track_data`}{Coordinates, `time_mins`, `n_locs`, timestamps,
#'       the computed metric `crs`, the caller's `orig_crs` (so the UD
#'       layer can return the UD in it), and the resolved per-fix
#'       `location_error`.}
#'   }
#'
#'   For multi-track input: a named list of `mt_dbgb_variance` objects.
#'   Tracks with too few locations are skipped with a warning.
#'
#' @details
#' The method extends the dBBMM by decomposing the variance into a
#' parallel component (along the direction of travel between consecutive
#' locations) and an orthogonal component (perpendicular to it). This
#' allows distinguishing directed movement (high parallel, low orthogonal
#' variance) from random movement (similar variances in both directions).
#'
#' A directionality index can be computed from the variances:
#' `I_d = (para - orth) / (para + orth)`, where values near 0 indicate
#' Brownian motion and positive values indicate directional movement.
#'
#' @references
#' Kranstauber, B., Safi, K., & Bartumeus, F. (2014). Bivariate Gaussian
#' bridges: directional factorization of diffusion in Brownian bridge
#' models. *Movement Ecology*, 2(1), 5. \doi{10.1186/2051-3933-2-5}
#'
#' @seealso [mt_dbgb_ud()] to compute the UD, [mt_motion_variance()] to
#'   extract variances as a data.frame, [mt_dbbmm_variance()] for the
#'   isotropic variant.
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
#' var_obj <- mt_dbgb_variance(f1_proj, location_error = 25,
#'                              margin = 15, window_size = 31)
#' var_obj
#'
#' # Directionality index
#' mv <- mt_motion_variance(var_obj)
#' I_d <- (mv$para - mv$orth) / (mv$para + mv$orth)
#' plot(mt_time(f1_proj), I_d, type = "l",
#'      xlab = "Time", ylab = "Directionality index")
#' }
#'
#' @export
mt_dbgb_variance <- function(object, location_error, window_size, margin,
                              parallel = FALSE, cores = NULL,
                              location_error_na = "median") {
  ## Record the caller's CRS *before* .ensure_projected so the UD layer can
  ## return the UD in it (the CRS-stratification contract).
  orig_crs <- if (mt_is_move2(object)) sf::st_crs(object) else NULL
  object <- .ensure_projected(object, "dynamic directional Brownian-bridge variance estimation")
  .validate_move2(object)

  if (mt_n_tracks(object) > 1) {
    return(.dbgb_variance_multi(object, location_error, window_size, margin,
                                 location_error_na = location_error_na,
                                 parallel = parallel, cores = cores,
                                 orig_crs = orig_crs))
  }

  .dbgb_variance_single(object, location_error, window_size, margin,
                          location_error_na = location_error_na,
                          parallel = parallel, cores = cores,
                          orig_crs = orig_crs)
}

#' @keywords internal
.dbgb_variance_single <- function(object, location_error, window_size, margin,
                                   location_error_na = "median",
                                   parallel = FALSE, cores = NULL,
                                   orig_crs = NULL) {
  td <- .extract_track_data(object)
  n <- td$n_locs

  location_error <- .resolve_loc_err_for_variance(
    location_error, object, n,
    na_replace = location_error_na
  )
  if (is.null(location_error)) location_error <- rep(0, n)

  if (n < window_size) {
    rlang::abort(
      sprintf("window_size (%d) cannot be larger than the number of locations (%d).",
              window_size, n),
      class = "move2utils_dbgb_variance_dyn_window_too_large")
  }
  if (any((c(margin, window_size) %% 2) != 1)) {
    rlang::abort("margin and window_size must both be odd.",
                 class = "move2utils_dbgb_variance_dyn_bad_window_parity")
  }

  ## Single-call C path: the per-window 4-step search + RMS
  ## aggregation runs entirely in C with optional OpenMP
  ## parallelisation over windows (see src/bgb_var_window_c.c ->
  ## bgb_var_break_track_c).  Numerics are identical to the per-
  ## window kernel (Brent's method per axis unchanged); the R-level
  ## lapply + do.call(rbind, ...) + aggregate() pipeline is replaced
  ## with one .Call.  `parallel` / `cores` are kept for API
  ## compatibility but are now ignored at this level.
  res <- .Call("bgb_var_break_track_c",
               td$x, td$y, td$time_mins, location_error,
               as.integer(window_size), as.integer(margin))
  para_sd <- res$para_sd
  orth_sd <- res$orth_sd
  n_estim <- ifelse(res$n_estim > 0L, as.numeric(res$n_estim), NA_real_)

  seg_interest <- !is.na(n_estim) & (n_estim == max(n_estim, na.rm = TRUE))

  structure(
    list(
      para_sd = para_sd,
      orth_sd = orth_sd,
      n_estim = n_estim,
      seg_interest = seg_interest,
      margin = margin,
      window_size = window_size,
      track_data = list(
        x = td$x, y = td$y,
        time_mins = td$time_mins,
        n_locs = n,
        crs = st_crs(object),
        orig_crs = if (is.null(orig_crs)) st_crs(object) else orig_crs,
        timestamps = mt_time(object),
        location_error = location_error
      )
    ),
    class = "mt_dbgb_variance"
  )
}

#' @keywords internal
.dbgb_variance_multi <- function(object, location_error, window_size, margin,
                                  location_error_na = "median",
                                  parallel = FALSE, cores = NULL,
                                  orig_crs = NULL) {
  ## Resolve column-name / "auto" inputs once on the whole multi-track
  ## object so the per-track slice is a clean numeric vector.  Scalar
  ## and NULL pass through unchanged.  Numeric vectors are
  ## length-checked against the full object.
  ids_all <- as.character(mt_track_id(object))
  if (!is.null(location_error) &&
      (is.character(location_error) ||
       (is.numeric(location_error) && length(location_error) > 1L))) {
    location_error <- .resolve_location_error(location_error, object,
                                               n = length(ids_all))
  }

  tracks <- .split_tracks(object)
  results <- list()
  skipped <- character(0)

  for (nm in names(tracks)) {
    trk <- tracks[[nm]]
    if (nrow(trk) < window_size) {
      skipped <- c(skipped, nm)
      next
    }
    le_trk <- if (is.numeric(location_error) && length(location_error) > 1L) {
      location_error[ids_all == nm]
    } else {
      location_error
    }
    results[[nm]] <- .dbgb_variance_single(
      trk, location_error = le_trk,
      margin = margin, window_size = window_size,
      location_error_na = location_error_na,
      parallel = parallel, cores = cores,
      orig_crs = orig_crs
    )
  }

  if (length(skipped) > 0) {
    rlang::warn(
      sprintf("Tracks skipped (fewer locations than window_size): %s",
              paste(skipped, collapse = ", ")),
      class = "move2utils_dbgb_variance_dyn_tracks_skipped")
  }
  if (length(results) == 0) {
    rlang::abort("No tracks had enough locations for the given window_size.",
                 class = "move2utils_dbgb_variance_dyn_all_too_short")
  }

  results
}

#' Print method for mt_dbgb_variance
#'
#' @param x An `mt_dbgb_variance` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mt_dbgb_variance <- function(x, ...) {
  cat("Dynamic Bivariate Gaussian Bridge \u2014 variance estimate\n")
  cat(sprintf("  Locations: %d\n", x$track_data$n_locs))
  cat(sprintf("  Window size: %d, Margin: %d\n", x$window_size, x$margin))
  cat(sprintf("  Parallel SD range: %.2f \u2013 %.2f\n",
              min(x$para_sd, na.rm = TRUE), max(x$para_sd, na.rm = TRUE)))
  cat(sprintf("  Orthogonal SD range: %.2f \u2013 %.2f\n",
              min(x$orth_sd, na.rm = TRUE), max(x$orth_sd, na.rm = TRUE)))
  invisible(x)
}

#' @rdname mt_motion_variance
#' @export
mt_motion_variance.mt_dbgb_variance <- function(x, ...) {
  data.frame(para = x$para_sd^2, orth = x$orth_sd^2)
}

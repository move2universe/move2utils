#' Dynamic Brownian Motion Variance Estimation
#'
#' Estimate the dynamic Brownian motion variance using a sliding window
#' approach with BIC-based breakpoint detection. Accepts single-track or
#' multi-track `move2` objects.
#'
#' @param object A `move2` object. A projected (metric) CRS is used as
#'   supplied; longitude/latitude input is auto-projected to a local AEQD
#'   internally (the variance values are CRS-invariant metres either way).
#'   Both single-track and multi-track objects are accepted. Empty
#'   geometries must be removed beforehand.
#' @param location_error Per-fix horizontal 1-sigma positional error,
#'   **in metres** (the AEQD / projected-CRS unit). The variance
#'   estimator treats this as the anchor measurement-error prior in the
#'   Brownian-bridge model.
#'
#'   Accepts:
#'   \itemize{
#'     \item \code{NULL} -- no per-fix error (treats fix positions as
#'           exact, equivalent to the legacy \code{location_error = 0}).
#'     \item a non-negative numeric scalar -- uniform sigma applied to
#'           every fix.
#'     \item a numeric vector of length \code{nrow(object)} -- per-fix
#'           sigma in metres.  Must match the row order of
#'           \code{object}; multi-track inputs are sliced per-track
#'           internally.
#'     \item a single character string -- name of a column in
#'           \code{object} containing per-fix sigma in metres
#'           (e.g. \code{"eobs_horizontal_accuracy_estimate"}).
#'     \item the literal \code{"auto"} -- probes
#'           \code{eobs_horizontal_accuracy_estimate} first, then
#'           \code{argos_lc} (mapped via the standard CLS / Vincent
#'           et al. 2002 class-to-sigma table). \code{gps_hdop} and
#'           similar quality indicators are deliberately not
#'           auto-resolved -- conversion to metres requires a
#'           sensor-specific multiplier the package cannot infer.
#'   }
#'   Negative values are rejected. \code{NA}s are imputed via
#'   \code{location_error_na} below.
#' @param location_error_na Strategy for filling \code{NA}s in a per-fix
#'   \code{location_error} vector (occurs when some fixes carry no
#'   quality info but others do).  One of \code{"median"} (default;
#'   per-track median of non-NA values -- conservative central tendency,
#'   robust), \code{"mean"} (per-track mean), \code{"zero"} (treat
#'   unknown errors as zero -- fakes certainty, only use when
#'   deliberately disregarding measurement error at those fixes), or
#'   \code{"approx"} (linear interpolation across positional index, with
#'   edge fill via nearest non-NA).
#' @param window_size Integer (must be odd). The number of locations in each
#'   sliding window. Larger values produce smoother variance estimates but
#'   miss short-term behavioural changes.
#' @param margin Integer (must be odd). The minimum number of locations on
#'   each side of a potential breakpoint within a window.
#'   Must satisfy `window_size >= 2 * margin + 1`.
#' @param parallel,cores Retained for API compatibility. Since 0.2.0
#'   the variance kernel is implemented in C with OpenMP threading
#'   (commit `bf8eb34`); these arguments are no-ops at the R level.
#'   To control thread count, set the `OMP_NUM_THREADS` environment
#'   variable before invoking R.
#'
#' @return For a single-track input: an `mt_dbbmm_variance` S3 object
#'   containing:
#'   \describe{
#'     \item{`variance`}{Numeric vector of estimated BM variances per
#'       location (`NA` for positions outside the estimable range at the
#'       track margins).}
#'     \item{`in_windows`}{Numeric vector: number of overlapping windows
#'       each location was estimated in.}
#'     \item{`interest`}{Logical vector: `TRUE` for locations fully
#'       covered by the sliding window (maximum overlap).}
#'     \item{`break_list`}{Integer vector of positions where behavioural
#'       breakpoints were detected.}
#'     \item{`window_size`, `margin`}{The parameters used.}
#'     \item{`track_data`}{List with `x`, `y`, `time_mins`, `n_locs`,
#'       `crs` (the metric CRS the variance was computed in), `orig_crs`
#'       (the caller's input CRS, used by the UD layer to return the UD in
#'       it), `timestamps`, and the resolved per-fix `location_error`.}
#'   }
#'
#'   For multi-track input: a named list of `mt_dbbmm_variance` objects,
#'   one per track. Tracks with fewer locations than `window_size` are
#'   skipped with a warning.
#'
#' @details
#' The method estimates Brownian motion variance using a leave-one-out
#' likelihood approach within a sliding window (Horne et al. 2007). At
#' each window position, it tests whether a single variance or two
#' variances (split at a breakpoint) better fit the data, using BIC for
#' model selection (Kranstauber et al. 2012). The window slides one
#' position at a time, and the final variance for each location is the
#' mean across all windows that covered it.
#'
#' The variance estimation and breakpoint testing are implemented in C
#' (Brent's method optimizer) for performance. The sliding window loop
#' can optionally be parallelised across CPU cores.
#'
#' @references
#' Horne, J. S., Garton, E. O., Krone, S. M., & Lewis, J. S. (2007).
#' Analyzing animal movements using Brownian bridges. *Ecology*, 88(9),
#' 2354-2363. \doi{10.1890/06-0957.1}
#'
#' Kranstauber, B., Kays, R., LaPoint, S. D., Wikelski, M., & Safi, K.
#' (2012). A dynamic Brownian bridge movement model to estimate utilization
#' distributions for heterogeneous animal movement. *Journal of Animal
#' Ecology*, 81(4), 738-746. \doi{10.1111/j.1365-2656.2012.01955.x}
#'
#' @seealso [mt_dbbmm_ud()] to compute the utilisation distribution from
#'   the variance estimate, [mt_motion_variance()] to extract the variance
#'   vector, [mt_dbgb_variance()] for the directional (bivariate) variant.
#'
#' @examples
#' \dontrun{
#' library(move2)
#' library(sf)
#'
#' # Load and prepare example data
#' fishers <- mt_read(mt_example())
#' fishers <- fishers[!st_is_empty(fishers), ]
#'
#' # Single track
#' f1 <- fishers[mt_track_id(fishers) == "F1", ]
#' f1_proj <- st_transform(f1, mt_aeqd_crs(f1))
#' var_obj <- mt_dbbmm_variance(f1_proj, location_error = 25,
#'                               window_size = 31, margin = 11)
#' var_obj
#' plot(mt_time(f1_proj), mt_motion_variance(var_obj),
#'      type = "l", xlab = "Time", ylab = "BM variance")
#'
#' # Multiple tracks
#' fishers_proj <- st_transform(fishers, mt_aeqd_crs(fishers))
#' var_list <- mt_dbbmm_variance(fishers_proj, location_error = 25,
#'                                window_size = 31, margin = 11)
#' names(var_list)
#' }
#'
#' @export
mt_dbbmm_variance <- function(object, location_error, window_size, margin,
                               parallel = FALSE, cores = NULL,
                               location_error_na = "median") {
  ## Record the caller's CRS *before* .ensure_projected so the UD layer can
  ## return the UD in it (the CRS-stratification contract).
  orig_crs <- if (mt_is_move2(object)) sf::st_crs(object) else NULL
  object <- .ensure_projected(object, "dynamic Brownian-bridge variance estimation")
  .validate_move2(object)

  if (mt_n_tracks(object) > 1) {
    return(.dbbmm_variance_multi(object, location_error, window_size, margin,
                                  location_error_na = location_error_na,
                                  parallel = parallel, cores = cores,
                                  orig_crs = orig_crs))
  }

  .dbbmm_variance_single(object, location_error, window_size, margin,
                           location_error_na = location_error_na,
                           parallel = parallel, cores = cores,
                           orig_crs = orig_crs)
}

#' @keywords internal
.dbbmm_variance_single <- function(object, location_error, window_size, margin,
                                    location_error_na = "median",
                                    parallel = FALSE, cores = NULL,
                                    orig_crs = NULL) {
  td <- .extract_track_data(object)
  n <- td$n_locs

  time_lag <- c(diff(td$time_mins), 0)
  location_error <- .resolve_loc_err_for_variance(
    location_error, object, n,
    na_replace = location_error_na
  )
  if (is.null(location_error)) location_error <- rep(0, n)

  if (n < window_size) {
    rlang::abort(
      sprintf("window_size (%d) cannot be larger than the number of locations (%d).",
              window_size, n),
      class = "move2utils_dbbmm_variance_dyn_window_too_large")
  }
  if (any((c(margin, window_size) %% 2) != 1)) {
    rlang::abort("margin and window_size must both be odd.",
                 class = "move2utils_dbbmm_variance_dyn_bad_window_parity")
  }
  if (window_size < 2 * margin + 1) {
    rlang::abort("window_size must be at least 2 * margin + 1.",
                 class = "move2utils_dbbmm_variance_dyn_window_too_small")
  }

  ## Single-call C path: the per-window sliding sweep + aggregation runs
  ## entirely in C with optional OpenMP parallelisation over windows
  ## (see src/bm_variance_c.c -> bm_variance_track_c).  The legacy
  ## R-level loop with .Call("bm_variance_window_c") per window plus
  ## do.call(rbind, ...) + aggregate() is replaced; numerics are
  ## identical (Brent's method per window unchanged), runtime drops
  ## by ~25-40% from the eliminated R glue plus near-linear scaling
  ## from OpenMP.  The `parallel` / `cores` arguments are retained for
  ## API compatibility but are now ignored at this level.
  res <- .Call("bm_variance_track_c",
               td$x, td$y, time_lag, location_error,
               as.integer(window_size), as.integer(margin))
  variance     <- res$variance
  in_windows   <- ifelse(res$n_estim > 0L, as.numeric(res$n_estim), NA_real_)
  breaks_found <- res$break_pos

  interest <- rep(FALSE, n)
  if (any(res$n_estim > 0L)) {
    interest[res$n_estim == max(res$n_estim)] <- TRUE
  }

  if (length(breaks_found) == 0L) breaks_found <- integer(0)

  structure(
    list(
      variance = variance,
      in_windows = in_windows,
      interest = interest,
      break_list = breaks_found,
      window_size = window_size,
      margin = margin,
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
    class = "mt_dbbmm_variance"
  )
}

#' @keywords internal
.dbbmm_variance_multi <- function(object, location_error, window_size, margin,
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
    n <- nrow(trk)
    if (n < window_size) {
      skipped <- c(skipped, nm)
      next
    }
    le_trk <- if (is.numeric(location_error) && length(location_error) > 1L) {
      location_error[ids_all == nm]
    } else {
      location_error
    }
    results[[nm]] <- .dbbmm_variance_single(
      trk, location_error = le_trk,
      window_size = window_size, margin = margin,
      location_error_na = location_error_na,
      parallel = parallel, cores = cores,
      orig_crs = orig_crs
    )
  }

  if (length(skipped) > 0) {
    rlang::warn(
      sprintf("Tracks skipped (fewer locations than window_size): %s",
              paste(skipped, collapse = ", ")),
      class = "move2utils_dbbmm_variance_dyn_tracks_skipped")
  }
  if (length(results) == 0) {
    rlang::abort("No tracks had enough locations for the given window_size.",
                 class = "move2utils_dbbmm_variance_dyn_all_too_short")
  }

  results
}

#' @keywords internal
.run_lapply <- function(X, FUN, parallel, cores) {
  if (parallel) {
    if (.Platform$OS.type == "windows") {
      message("Note: parallel processing uses mclapply which is not available ",
              "on Windows. Falling back to sequential processing.")
      lapply(X, FUN)
    } else {
      if (is.null(cores)) cores <- max(1, parallel::detectCores() - 1)
      parallel::mclapply(X, FUN, mc.cores = cores)
    }
  } else {
    lapply(X, FUN)
  }
}


#' Print method for mt_dbbmm_variance
#'
#' @param x An `mt_dbbmm_variance` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mt_dbbmm_variance <- function(x, ...) {
  cat("Dynamic Brownian Bridge Movement Model \u2014 variance estimate\n")
  cat(sprintf("  Locations: %d\n", x$track_data$n_locs))
  cat(sprintf("  Window size: %d, Margin: %d\n", x$window_size, x$margin))
  cat(sprintf("  Breakpoints found: %d\n", length(x$break_list)))
  cat(sprintf("  Variance range: %.2f \u2013 %.2f\n",
              min(x$variance, na.rm = TRUE), max(x$variance, na.rm = TRUE)))
  invisible(x)
}


#' Extract Motion Variance
#'
#' Extract the estimated movement variance from a dBBMM or dBGB variance
#' object.
#'
#' @param x An `mt_dbbmm_variance` or `mt_dbgb_variance` object, as
#'   returned by [mt_dbbmm_variance()] or [mt_dbgb_variance()].
#' @param ... Ignored.
#'
#' @return For `mt_dbbmm_variance`: a numeric vector of variances (one
#'   per location, `NA` at track margins).
#'
#'   For `mt_dbgb_variance`: a `data.frame` with columns `para`
#'   (parallel variance) and `orth` (orthogonal variance).
#'
#' @seealso [mt_dbbmm_variance()], [mt_dbgb_variance()]
#'
#' @examples
#' \dontrun{
#' library(move2)
#' library(sf)
#' fishers <- mt_read(mt_example())
#' fishers <- fishers[!st_is_empty(fishers), ]
#' f1 <- fishers[mt_track_id(fishers) == "F1", ]
#' f1_proj <- st_transform(f1, mt_aeqd_crs(f1))
#' var_obj <- mt_dbbmm_variance(f1_proj, location_error = 25,
#'                               window_size = 31, margin = 11)
#' head(mt_motion_variance(var_obj))
#' }
#'
#' @export
mt_motion_variance <- function(x, ...) {
  UseMethod("mt_motion_variance")
}

#' @rdname mt_motion_variance
#' @export
mt_motion_variance.mt_dbbmm_variance <- function(x, ...) {
  x$variance
}

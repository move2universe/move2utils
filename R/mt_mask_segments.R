#' Exclude segments from a dBBMM or dBGB utilisation distribution
#'
#' Mark specific track segments so that their Brownian bridges are
#' excluded when [mt_dbbmm_ud()] or [mt_dbgb_ud()] computes the
#' utilisation distribution. A segment `i` is the bridge from location
#' `i` to location `i + 1`; setting it to excluded drops the
#' contribution of that bridge while leaving the variance estimates and
#' the neighbouring bridges untouched.
#'
#' The helper is a thin dispatcher over the per-class include-masks
#' (`interest` on `mt_dbbmm_variance`, `seg_interest` on
#' `mt_dbgb_variance`). Both masks have identical semantics: `TRUE`
#' means "include this bridge in the UD", `FALSE` means "exclude it".
#' Using this helper keeps user code identical across the two variance
#' classes.
#'
#' The common use case is removing long-gap bridges that would otherwise
#' fan the UD out across regions the animal never visited. See
#' `vignette("UD_gap_aware_ud", package = "move2utils")` for a worked
#' example on the fisher data.
#'
#' @param x A variance object of class `mt_dbbmm_variance` or
#'   `mt_dbgb_variance`, or a named list of such objects as returned by
#'   multi-track dispatch.
#' @param segments Integer vector of segment indices to exclude (each in
#'   `1:(n_locs - 1)`), **or** a logical vector of length `n_locs - 1`
#'   or `n_locs` where `TRUE` means "exclude". `NA` entries are ignored.
#'
#' @return The variance object (or list of objects) with the mask field
#'   updated.
#'
#' @examples
#' \dontrun{
#' library(move2)
#' library(sf)
#' library(units)
#'
#' fishers <- mt_read(mt_example())
#' fishers <- fishers[!st_is_empty(fishers), ]
#' leroy   <- fishers[mt_track_id(fishers) == "M1", ]
#' leroy   <- st_transform(leroy, mt_aeqd_crs(leroy))
#'
#' lag_h <- as.numeric(set_units(mt_time_lags(leroy), "hour"))
#' long_gap <- which(lag_h > 6)
#'
#' var <- mt_dbbmm_variance(leroy, location_error = 25,
#'                           window_size = 31, margin = 11)
#' var <- mt_mask_segments(var, long_gap)
#' ud  <- mt_dbbmm_ud(var, location_error = 25, raster = 100, ext = 0.5)
#' }
#'
#' @seealso [mt_dbbmm_variance()], [mt_dbgb_variance()],
#'   [mt_dbbmm_ud()], [mt_dbgb_ud()].
#' @export
mt_mask_segments <- function(x, segments) UseMethod("mt_mask_segments")

#' @export
mt_mask_segments.mt_dbbmm_variance <- function(x, segments) {
  .apply_segment_mask(x, segments, field = "interest")
}

#' @export
mt_mask_segments.mt_dbgb_variance <- function(x, segments) {
  .apply_segment_mask(x, segments, field = "seg_interest")
}

#' @export
mt_mask_segments.list <- function(x, segments) {
  if (is.list(segments) && !is.null(names(segments)) &&
      all(names(segments) %in% names(x))) {
    out <- x
    for (nm in names(segments)) {
      out[[nm]] <- mt_mask_segments(out[[nm]], segments[[nm]])
    }
    return(out)
  }
  lapply(x, mt_mask_segments, segments = segments)
}

#' @export
mt_mask_segments.default <- function(x, segments) {
  rlang::abort(
    sprintf(paste0(
      "`mt_mask_segments()` is defined for objects of class ",
      "`mt_dbbmm_variance` or `mt_dbgb_variance` (or a list of these). ",
      "Got class: %s."),
      paste(class(x), collapse = "/")),
    class = "move2utils_mt_mask_segments_unsupported_class")
}

#' @keywords internal
.apply_segment_mask <- function(x, segments, field) {
  n <- length(x[[field]])

  if (is.logical(segments)) {
    if (!length(segments) %in% c(n - 1L, n)) {
      rlang::abort(
        sprintf("logical `segments` must have length n-1 (%d) or n (%d); got %d.",
                n - 1L, n, length(segments)),
        class = "move2utils_mt_mask_segments_bad_length")
    }
    segments[is.na(segments)] <- FALSE
    segments <- which(segments)
  } else if (is.numeric(segments)) {
    segments <- segments[!is.na(segments)]
    if (length(segments) && (any(segments < 1) || any(segments > n - 1L))) {
      rlang::abort(
        sprintf("`segments` must be integers in 1:%d.", n - 1L),
        class = "move2utils_mt_mask_segments_out_of_range")
    }
    segments <- as.integer(segments)
  } else {
    rlang::abort("`segments` must be an integer or logical vector.",
                 class = "move2utils_mt_mask_segments_bad_type")
  }

  if (length(segments)) {
    x[[field]][segments] <- FALSE
  }
  x
}

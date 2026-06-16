#' Filter fixes by GPS-geometry quality
#'
#' Pre-processing step that removes locations whose GPS geometry is
#' unreliable before outlier-probability scoring. Bad-geometry fixes
#' (too few satellites, high DOP, large eObs horizontal-accuracy
#' estimate) typically land far from the true position and can dominate
#' the histogram range used by \code{\link{mt_flag_outliers}}.
#' Satellite count is the most canonical of these indicators: GPS
#' requires at least four satellites for a trilateration-based fix, and
#' fixes at the four-satellite boundary have pathological error
#' geometry. A threshold of \code{sat_min = 5} matches standard
#' telemetry practice.
#'
#' Columns are detected by Movebank conventions
#' (\code{gps_satellite_count}, \code{gps_dop},
#' \code{eobs_horizontal_accuracy_estimate}). Missing columns are
#' skipped; to disable a criterion on a track that has the column,
#' pass \code{NULL} for the corresponding argument.
#'
#' Empty geometries (rows whose location is missing) are the limit
#' case of an untrustworthy fix and are dropped by default. They
#' have no bearing on sat/DOP/hacc but cannot be scored by any
#' downstream detector, and their presence contaminates neighbour-
#' based analyses. Set \code{drop_empty = FALSE} to keep them.
#'
#' @param x A \code{move2} object.
#' @param sat_min Integer. Keep fixes with at least this many
#'   satellites. Default 5. Set to \code{NULL} to skip.
#' @param dop_max Numeric. Keep fixes with DOP at or below this value.
#'   Default 10. Set to \code{NULL} to skip.
#' @param hacc_max Numeric. Keep fixes whose
#'   \code{eobs_horizontal_accuracy_estimate} (metres) is at or below
#'   this value. Default 100. Set to \code{NULL} to skip.
#' @param drop_empty Logical. If \code{TRUE} (default), drop fixes
#'   whose geometry is empty (missing location). Set to \code{FALSE}
#'   to keep them -- in which case downstream detectors will handle
#'   them as non-scorable rows.
#' @param verbose Logical. If \code{TRUE} (default), report per-criterion
#'   drop counts via \code{message()}.
#'
#' @return A \code{move2} object with unreliable fixes removed. A
#'   message summarises how many were dropped by each criterion.
#'
#' @examples
#' \dontrun{
#' library(move2)
#' ## a GPS track downloaded from Movebank, carrying satellite-count /
#' ## DOP / horizontal-accuracy columns
#' x <- movebank_download_study("Your study name", ...)
#' x <- x[!sf::st_is_empty(x), ]
#' x <- mt_filter_unique(x, criterion = "first")
#' x <- mt_filter_gps_quality(x)
#' }
#'
#' @export
mt_filter_gps_quality <- function(x, sat_min = 5, dop_max = 10,
                                   hacc_max = 100, drop_empty = TRUE,
                                   verbose = TRUE) {

  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }

  n0 <- nrow(x)
  keep <- rep(TRUE, n0)
  report <- list()
  any_quality_col <- FALSE

  ## empty geometries: treated as the limit case of bad geometry.
  ## Dropped by default; can be disabled for callers that want to
  ## preserve row parity with the original object.
  if (isTRUE(drop_empty)) {
    empty <- sf::st_is_empty(x)
    drop  <- !is.na(empty) & empty
    keep  <- keep & !drop
    report[["empty"]] <- list(n = sum(drop),
                              label = "empty geometry")
  }

  ## Accept both name conventions: movebank_download_study returns
  ## underscore names; mt_read() on a raw Movebank CSV keeps dashes.
  find_col <- function(underscore_name) {
    if (underscore_name %in% names(x)) return(underscore_name)
    dashed <- gsub("_", "-", underscore_name, fixed = TRUE)
    if (dashed %in% names(x)) return(dashed)
    NA_character_
  }

  ## satellite count
  if (!is.null(sat_min)) {
    col <- find_col("gps_satellite_count")
    if (!is.na(col)) {
      any_quality_col <- TRUE
      v <- suppressWarnings(as.numeric(x[[col]]))
      drop <- !is.na(v) & v < sat_min
      keep <- keep & !drop
      report[["sat"]] <- list(n = sum(drop),
                              label = sprintf("sat < %d", as.integer(sat_min)))
    }
  }

  ## dilution of precision
  if (!is.null(dop_max)) {
    col <- find_col("gps_dop")
    if (!is.na(col)) {
      any_quality_col <- TRUE
      v <- suppressWarnings(as.numeric(x[[col]]))
      drop <- !is.na(v) & v > dop_max
      keep <- keep & !drop
      report[["dop"]] <- list(n = sum(drop),
                              label = sprintf("DOP > %g", dop_max))
    }
  }

  ## eObs horizontal-accuracy estimate (metres)
  if (!is.null(hacc_max)) {
    col <- find_col("eobs_horizontal_accuracy_estimate")
    if (!is.na(col)) {
      any_quality_col <- TRUE
      v <- suppressWarnings(as.numeric(x[[col]]))
      drop <- !is.na(v) & v > hacc_max
      keep <- keep & !drop
      report[["hacc"]] <- list(n = sum(drop),
                               label = sprintf("hacc > %g m", hacc_max))
    }
  }

  if (!any_quality_col && (!isTRUE(drop_empty) || report[["empty"]]$n == 0)) {
    if (verbose) {
      message("No GPS-quality columns found (expected one of ",
              "gps_satellite_count, gps_dop, ",
              "eobs_horizontal_accuracy_estimate). Returning unchanged.")
    }
    return(x)
  }

  n_kept <- sum(keep)
  n_dropped <- n0 - n_kept

  if (verbose) {
    message(sprintf(
      "GPS-quality filter: kept %d of %d fixes (%d dropped).",
      n_kept, n0, n_dropped))
    for (r in report) {
      if (r$n > 0) {
        message(sprintf("  %s: %d", r$label, r$n))
      }
    }
  }

  x[keep, ]
}

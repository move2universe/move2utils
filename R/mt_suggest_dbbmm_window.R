## Diagnostic helper for dBBMM / dBGB variance estimation.
##
## mt_dbbmm_variance() and mt_dbgb_variance() require two integer
## hyperparameters: window_size (number of fixes per sliding window)
## and margin (minimum fixes either side of any tested breakpoint
## within a window).  The defaults in the literature (window_size=31,
## margin=11) are reasonable for half-hourly mammal tracks, but on
## coarser or finer sampling the temporal coverage of the window
## can become inappropriate -- too short to bridge a behavioural
## episode, or too long to resolve one.
##
## This helper inspects the track's temporal sampling and suggests
## a window_size + margin pair whose temporal coverage matches a
## user-specified target (default 4 hours per window).  It does not
## decide for the user; the suggestion is one starting point among
## many.


#' Suggest dBBMM / dBGB sliding-window parameters from track sampling
#'
#' Inspects a \code{move2} object's temporal sampling and proposes a
#' \code{window_size} and \code{margin} pair for
#' \code{\link{mt_dbbmm_variance}} or \code{\link{mt_dbgb_variance}}
#' whose one-window temporal coverage matches a user-specified target.
#' Diagnostic only -- the function does not call the variance
#' estimator, just prints the suggestion (and optionally plots the
#' time-lag distribution).
#'
#' @details
#' The dBBMM / dBGB variance estimator slides a window of
#' \code{window_size} fixes along the track and tests within each
#' window whether a single Brownian-motion variance or two variances
#' (split at a breakpoint) better fit the data, requiring at least
#' \code{margin} fixes on each side of any tested breakpoint
#' (Kranstauber et al. 2012).  The temporal coverage of one window is
#' \code{window_size * median(dt)} for typical tracks.
#'
#' This helper picks \code{window_size} as the nearest odd integer to
#' \code{target_hours * 3600 / median(dt)}, clamped to a minimum of 11
#' (the smallest window that gives reasonable BIC stability) and to
#' at most \code{floor(n/4)} so several non-overlapping windows fit
#' along the track.  \code{margin} is set to the nearest odd integer
#' to \code{(window_size - 1) / 4}, satisfying the
#' \code{window_size >= 2 * margin + 1} constraint with comfortable
#' headroom.
#'
#' Heuristics aside, the right \code{window_size} depends on the
#' biology of the species and the question being asked.  Sweep a few
#' values around the suggestion when the result matters.
#'
#' @param x A \code{move2} object.  Single- or multi-track.  The CRS
#'   does not matter -- only time lags are used.
#' @param target_hours Numeric.  Target temporal coverage of one
#'   sliding window in hours.  Default 4.  The literature default of
#'   \code{window_size = 31, margin = 11} corresponds to ~4 h on
#'   half-hourly tracks and ~30 min on minute-spaced tracks.
#' @param plot Logical.  If \code{TRUE} (default), plot the time-lag
#'   distribution and overlay the suggested window's temporal span.
#'
#' @return Invisibly, a list with components
#'   \describe{
#'     \item{\code{window_size}}{Suggested integer window size.}
#'     \item{\code{margin}}{Suggested integer margin.}
#'     \item{\code{median_dt_secs}}{Median time lag in seconds.}
#'     \item{\code{target_hours}}{Target coverage used.}
#'     \item{\code{window_hours}}{Actual coverage of the suggested
#'           window in hours
#'           (\code{window_size * median_dt_secs / 3600}).}
#'     \item{\code{n}}{Number of locations.}
#'     \item{\code{n_windows}}{Approximate number of non-overlapping
#'           windows that fit in the track.}
#'   }
#'   For multi-track input the list is named per individual.
#'
#' @examples
#' \dontrun{
#' library(move2)
#' fishers <- mt_read(mt_example())
#' s <- mt_suggest_dbbmm_window(fishers, target_hours = 4)
#' s
#' ## use the suggestion:
#' var_obj <- mt_dbbmm_variance(fishers, location_error = 25,
#'                                window_size = s$window_size,
#'                                margin      = s$margin)
#' }
#'
#' @references
#' Kranstauber, B., Kays, R., LaPoint, S. D., Wikelski, M., & Safi, K.
#' (2012). A dynamic Brownian bridge movement model to estimate
#' utilization distributions for heterogeneous animal movement.
#' \emph{Journal of Animal Ecology}, 81(4), 738-746.
#' \doi{10.1111/j.1365-2656.2012.01955.x}
#'
#' @seealso \code{\link{mt_dbbmm_variance}},
#'   \code{\link{mt_dbgb_variance}}.
#'
#' @importFrom move2 mt_time_lags mt_track_id mt_n_tracks
#' @importFrom graphics par hist abline legend mtext
#' @export
mt_suggest_dbbmm_window <- function(x, target_hours = 4, plot = TRUE) {
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  if (!is.numeric(target_hours) || length(target_hours) != 1L ||
      is.na(target_hours) || target_hours <= 0) {
    rlang::abort("`target_hours` must be a positive scalar.",
                 class = "move2utils_mt_suggest_dbbmm_window_bad_target_hours")
  }

  ## ---- multi-individual dispatch ----
  ids <- move2::mt_track_id(x)
  unique_ids <- unique(ids)
  if (length(unique_ids) > 1L) {
    out <- lapply(unique_ids, function(id) {
      message(sprintf("--- %s ---", id))
      mt_suggest_dbbmm_window(x[ids == id, ],
                              target_hours = target_hours,
                              plot = FALSE)
    })
    names(out) <- unique_ids
    if (plot) .plot_dbbmm_window_suggestion_multi(x, out, target_hours)
    return(invisible(out))
  }

  n <- nrow(x)
  if (n < 11L) {
    message("Too few locations for dBBMM (<11). Returning NULL.")
    return(invisible(NULL))
  }

  dt_s <- as.numeric(move2::mt_time_lags(x, units = "secs"))
  dt_s <- dt_s[is.finite(dt_s) & dt_s > 0]
  if (length(dt_s) < 2L) {
    message("Too few finite time lags. Returning NULL.")
    return(invisible(NULL))
  }
  med_dt <- stats::median(dt_s)

  ## ---- suggest window_size ---------------------------------------
  ## Target: window covers `target_hours` of real time at the typical
  ## sampling cadence.  Round to odd; clamp to [11, floor(n/4) | odd].
  raw_w  <- target_hours * 3600 / med_dt
  cap_w  <- max(11L, .floor_odd(n / 4L))
  w_int  <- .nearest_odd(raw_w)
  w      <- max(11L, min(w_int, cap_w))

  ## ---- suggest margin --------------------------------------------
  ## (window - 1) / 4 rounded to nearest odd, with the
  ## 2*margin + 1 <= window constraint always satisfied.
  m_raw <- (w - 1L) / 4
  m_int <- max(3L, .nearest_odd(m_raw))
  m_max <- .floor_odd((w - 1L) / 2L)
  margin <- min(m_int, m_max)

  out <- list(
    window_size    = as.integer(w),
    margin         = as.integer(margin),
    median_dt_secs = unname(med_dt),
    target_hours   = target_hours,
    window_hours   = unname(w * med_dt / 3600),
    n              = n,
    n_windows      = floor(n / w)
  )
  class(out) <- "mt_dbbmm_window_suggestion"

  ## ---- report ----------------------------------------------------
  message(sprintf(
    "  n = %d, median dt = %s; suggested window_size = %d (covers ~%.2f h), margin = %d.",
    out$n,
    .format_dt(out$median_dt_secs),
    out$window_size, out$window_hours, out$margin))
  if (out$window_size == 11L && out$window_hours < target_hours / 2) {
    message("  Note: hit the window_size = 11 floor; track sampling is too coarse ",
            "for the target window coverage.")
  }
  if (out$n_windows < 4L) {
    message(sprintf(
      "  Note: only ~%d non-overlapping window(s) fit in the track; consider a smaller `target_hours`.",
      out$n_windows))
  }

  if (plot) .plot_dbbmm_window_suggestion(dt_s, out, target_hours)

  invisible(out)
}


#' @export
print.mt_dbbmm_window_suggestion <- function(x, ...) {
  cat("dBBMM / dBGB sliding-window suggestion\n")
  cat(sprintf("  window_size : %d (covers %.2f h)\n",
              x$window_size, x$window_hours))
  cat(sprintf("  margin      : %d\n", x$margin))
  cat(sprintf("  n locations : %d\n", x$n))
  cat(sprintf("  median dt   : %s\n", .format_dt(x$median_dt_secs)))
  cat(sprintf("  target      : %.2f h per window\n", x$target_hours))
  cat(sprintf("  approx %d non-overlapping windows fit in this track\n",
              x$n_windows))
  invisible(x)
}


## Round to nearest odd integer (>= 1).
##
## @keywords internal
.nearest_odd <- function(z) {
  k <- as.integer(round(z))
  if (k %% 2L == 0L) k <- k + 1L
  if (k < 1L) k <- 1L
  k
}

## Largest odd integer <= z (>= 1).
##
## @keywords internal
.floor_odd <- function(z) {
  k <- as.integer(floor(z))
  if (k %% 2L == 0L) k <- k - 1L
  if (k < 1L) k <- 1L
  k
}

## Pretty-print a duration in seconds.
##
## @keywords internal
.format_dt <- function(s) {
  if (s < 60)        return(sprintf("%.1f s",   s))
  if (s < 3600)      return(sprintf("%.1f min", s / 60))
  if (s < 86400)     return(sprintf("%.2f h",   s / 3600))
  sprintf("%.2f d", s / 86400)
}


## Diagnostic plot for a single-track suggestion.
##
## @keywords internal
.plot_dbbmm_window_suggestion <- function(dt_s, sugg, target_hours) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(4, 4, 3, 1))

  graphics::hist(log10(dt_s),
                 breaks = 30,
                 xlab = expression(log[10](Delta * t / s)),
                 main = "dBBMM window suggestion",
                 col = grDevices::adjustcolor("steelblue", 0.5),
                 border = "white")
  graphics::abline(v = log10(sugg$median_dt_secs),
                   lty = 2, col = "navy", lwd = 1.5)
  graphics::mtext(side = 3, line = 0.3, cex = 0.85,
                  sprintf("window_size = %d (~%.2f h), margin = %d",
                          sugg$window_size, sugg$window_hours, sugg$margin))
  graphics::legend("topright", bty = "n",
                   legend = c("median dt"),
                   lty = 2, col = "navy", lwd = 1.5)
  invisible(NULL)
}


## Multi-track plot: one panel per individual.
##
## @keywords internal
.plot_dbbmm_window_suggestion_multi <- function(x, suggs, target_hours) {
  ids <- names(suggs)
  k <- length(ids)
  if (k == 0L) return(invisible(NULL))
  ncol <- min(k, 3L)
  nrow <- ceiling(k / ncol)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1))
  for (id in ids) {
    s <- suggs[[id]]
    if (is.null(s)) {
      graphics::plot.new(); graphics::title(main = id)
      next
    }
    dt_s <- as.numeric(move2::mt_time_lags(
      x[move2::mt_track_id(x) == id, ], units = "secs"))
    dt_s <- dt_s[is.finite(dt_s) & dt_s > 0]
    if (length(dt_s) < 2L) {
      graphics::plot.new(); graphics::title(main = id)
      next
    }
    graphics::hist(log10(dt_s),
                   breaks = 25,
                   xlab = expression(log[10](Delta * t / s)),
                   main = id,
                   col = grDevices::adjustcolor("steelblue", 0.5),
                   border = "white")
    graphics::abline(v = log10(s$median_dt_secs),
                     lty = 2, col = "navy", lwd = 1.5)
    graphics::mtext(side = 3, line = 0.3, cex = 0.7,
                    sprintf("ws = %d, m = %d (%.2f h)",
                            s$window_size, s$margin, s$window_hours))
  }
  invisible(NULL)
}

#' Post-run diagnostic suite for mt_clean_track results
#'
#' Six diagnostic panels that let you eyeball whether an
#' \code{\link{mt_clean_track}} run was healthy.  Analogous to
#' \code{plot.lm()} for a linear model: each panel addresses a
#' distinct failure mode and prints a brief interpretive note
#' pointing at a remedy when something looks off.
#'
#' @details
#' Run \code{mt_diagnose_clean_track()} on the object returned by
#' \code{mt_clean_track(..., remove = FALSE)}.  The function reads the
#' flag columns and recomputes the small amount of additional
#' information needed (per-fix step speeds; nothing requires re-
#' running the per-fix detectors).  On multi-individual input the
#' first five panels focus on a single individual (default: the one
#' with the highest flag rate); a sixth panel summarises the cohort.
#'
#' \strong{The six panels}
#'
#' \describe{
#'   \item{1. log-speed density with detected modes}{KDE of
#'     \eqn{\log(\text{step speed})} with vertical lines for each
#'     substantive mode and the v_max cap that was used.  A single
#'     dominant mode + sparse upper tail = healthy.  Two or more
#'     substantive modes = bimodal behaviour (rest + flight); the
#'     per-fix detectors threshold against a single distribution and
#'     can over-flag the smaller mode.  Remedy: see
#'     \code{vignette('state_conditional')} (in preparation), or run
#'     each behavioural state separately.}
#'   \item{2. flag rate vs time}{Rolling-window flag rate over the
#'     timeline.  Roughly flat at <1\% with isolated spikes at known
#'     noise events = healthy.  A sustained elevated band over a
#'     contiguous window = the migration-over-flagging signature; the
#'     detectors are catching legitimate movement, not errors.
#'     Remedy: filter out that window or run it separately with
#'     stricter thresholds.}
#'   \item{3. per-detector activity history}{Bar chart of how many
#'     fixes each combination of detectors flagged during the
#'     iteration loop, BEFORE the conjunction rule was applied.
#'     "consensus" bins (>=2 detectors agree) are the high-confidence
#'     flags; "single-detector" bins are fixes one detector flagged
#'     that the conjunction subsequently rejected.  A large
#'     bridge-only or prob-only bar means that detector is firing
#'     noisily on the data; a track with mostly consensus bars is
#'     converging cleanly.}
#'   \item{4. cumulative flagging by iteration}{Cumulative flag count
#'     across the per-fix iteration loop.  Rapid plateau in 2-4
#'     iterations = healthy.  Linear growth without plateau = self-
#'     reinforcing flagging; the detectors are not converging on a
#'     stable flag set.}
#'   \item{5. consecutive-flag run length distribution}{Histogram of
#'     run lengths of consecutively flagged fixes.  Most flags
#'     isolated (length 1) = discrete errors, the case the package is
#'     calibrated for.  A heavy tail of long runs without
#'     \code{error_class = "block"} = a sustained behavioural state
#'     being mistaken for outliers.}
#'   \item{6. per-individual flag rates (multi-track only)}{Dot plot
#'     of flag rate per individual with reference lines at 0.5\%
#'     (typical clean), 2\% (suspect), 5\% (likely state-bimodal).
#'     Use this to triage which individuals need further attention.}
#' }
#'
#' @param x A \code{move2} object returned by
#'   \code{\link{mt_clean_track}} with \code{remove = FALSE}.  Must
#'   carry the flag columns the orchestrator attaches
#'   (\code{is_outlier}, \code{flagged_by_*}, \code{flag_iteration},
#'   \code{error_class}).
#' @param individual For multi-track input, the track id to focus the
#'   first five panels on.  Default \code{NULL} picks the individual
#'   with the highest flag rate (the one most worth diagnosing).
#'   Pass any track id to override.
#' @param window_days Numeric.  Target width of the rolling window (in
#'   days) used by Panel 2.  Default 7.  Larger windows smooth more;
#'   smaller windows resolve finer time-localised flag clusters.  On
#'   short tracks the width is automatically shrunk so the panel always
#'   has at least a handful of windows to plot (the previous fixed
#'   7-day window produced a "too few windows" placeholder on tracks
#'   spanning under two weeks); the width actually used is annotated on
#'   the panel.
#' @param cex_scale Numeric multiplier (default 1.2) applied to plot
#'   titles, axis labels, tick labels, legends, and in-panel
#'   annotations.  The six-panel \code{mfrow} layout shrinks base text
#'   substantially; bump this if the writing is hard to read after
#'   saving to a file, or set to 1 for the historical sizing.
#' @param plot Logical.  If \code{TRUE} (default), render the six-
#'   panel figure.  Set \code{FALSE} to skip plotting and just receive
#'   the diagnostic data.
#' @param silent Logical.  If \code{FALSE} (default), print the
#'   interpretive notes for each panel that tripped a concern.  Set
#'   \code{TRUE} to suppress.
#'
#' @return Invisibly, a list with components \code{by_individual}
#'   (data.frame of per-track flag rates), \code{run_lengths} (integer
#'   vector of consecutive-flag run lengths on the focused track),
#'   \code{modes} (numeric vector of detected substantive modes in
#'   m/s), and \code{notes} (character vector of interpretive
#'   messages emitted).  These let downstream code consume the
#'   diagnostic without re-running the function.
#'
#' @section What the diagnostic does NOT do:
#' It does not re-run the per-fix detectors and therefore cannot show
#' the underlying bridge \eqn{\eta} or joint-probability distributions
#' (a possible follow-up).  It also assumes a single behavioural-state
#' threshold was applied; tracks where the user has already segmented
#' by state and run \code{mt_clean_track} per segment will read
#' "healthy" on every panel because each segment IS unimodal.
#'
#' @examples
#' \dontrun{
#' library(move2)
#' x <- mt_read(system.file("extdata/synthetic_tracks.csv.gz",
#'                            package = "move2utils"))
#' x <- x[!sf::st_is_empty(x), ]
#' res <- mt_clean_track(x, plot = FALSE, remove = FALSE)
#' mt_diagnose_clean_track(res)
#' }
#'
#' @seealso \code{\link{mt_clean_track}} for the orchestrator that
#'   produces the input.
#'
#' @importFrom move2 mt_time mt_track_id mt_time_lags
#' @importFrom sf st_coordinates
#' @importFrom stats density
#' @importFrom graphics par plot points lines abline legend axis
#'   barplot mtext rect text title
#' @importFrom grDevices adjustcolor
#' @export
mt_diagnose_clean_track <- function(x,
                                     individual = NULL,
                                     window_days = 7,
                                     cex_scale = 1.2,
                                     plot = TRUE,
                                     silent = FALSE) {
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object returned by mt_clean_track().",
                 class = "move2utils_input_not_move2")
  }
  required_cols <- c("is_outlier", "flagged_by_bridge",
                     "flagged_by_prob", "flagged_by_speed",
                     "flag_iteration", "error_class")
  miss <- setdiff(required_cols, names(x))
  if (length(miss)) {
    rlang::abort(
      sprintf("`x` is missing flag columns from mt_clean_track(): %s.\n  Did you call mt_clean_track(..., remove = FALSE)?",
              paste(miss, collapse = ", ")),
      class = "move2utils_mt_diagnose_clean_track_missing_flag_columns")
  }

  say <- function(...) if (!silent) message(...)

  ids <- as.character(move2::mt_track_id(x))
  unique_ids <- unique(ids)
  is_multi <- length(unique_ids) > 1L

  ## --- per-individual flag rates (used by Panel 6 + focus picking) ---
  by_individual <- data.frame(
    individual = unique_ids,
    n          = vapply(unique_ids, function(id) sum(ids == id), integer(1)),
    flagged    = vapply(unique_ids, function(id)
                          sum(x$is_outlier[ids == id], na.rm = TRUE),
                          integer(1)),
    stringsAsFactors = FALSE
  )
  by_individual$pct <- 100 * by_individual$flagged / pmax(by_individual$n, 1L)
  by_individual <- by_individual[order(-by_individual$pct), ]

  ## --- pick the focused individual ---
  if (is.null(individual)) {
    individual <- by_individual$individual[1]
  }
  if (!individual %in% unique_ids) {
    rlang::abort(
      sprintf("`individual = '%s'` not found in track ids: %s",
              individual, paste(unique_ids, collapse = ", ")),
      class = "move2utils_mt_diagnose_clean_track_unknown_individual")
  }
  focus_idx <- ids == individual
  xf <- x[focus_idx, ]

  ## --- accumulators for the return value ---
  notes <- character()

  ## --- prepare panel data --------------------------------------------
  ## Panel 1: log-speed density with substantive modes
  step_m <- .step_lengths_fast(xf)
  dt_s   <- as.numeric(move2::mt_time_lags(xf, units = "secs"))
  speed  <- ifelse(is.finite(step_m) & is.finite(dt_s) & dt_s > 0,
                   step_m / dt_s, NA_real_)
  v_pos  <- speed[is.finite(speed) & speed > 0]
  v_max_used <- attr(x, "v_max_used")
  if (is.list(v_max_used) || length(v_max_used) > 1L) {
    v_max_focus <- as.numeric(v_max_used[as.character(individual)])
    if (length(v_max_focus) == 0L || is.na(v_max_focus)) v_max_focus <- NA_real_
  } else if (length(v_max_used) == 1L) {
    v_max_focus <- as.numeric(v_max_used)
  } else {
    v_max_focus <- NA_real_
  }

  modes <- .find_substantive_modes(v_pos)
  if (length(modes) >= 2L) {
    notes <- c(notes, sprintf(
      "Panel 1: %d substantive modes detected at %s m/s -- bimodal behaviour. The per-fix detectors threshold against a single distribution; consider state-conditional analysis or filtering to one mode before cleaning.",
      length(modes),
      paste(sprintf("%.2f", modes), collapse = ", ")))
  }
  if (is.finite(v_max_focus) && length(modes) > 0L &&
      v_max_focus < max(modes)) {
    notes <- c(notes, sprintf(
      "Panel 1: v_max (%.2f m/s) sits BELOW the rightmost substantive mode (%.2f m/s). The cap is severing a real activity mode, not isolating outliers.",
      v_max_focus, max(modes)))
  }

  ## Panel 2: flag rate vs time (rolling window)
  ts <- as.numeric(move2::mt_time(xf))   # seconds since epoch
  is_out <- xf$is_outlier
  panel2 <- .rolling_flag_rate(ts, is_out, window_secs = window_days * 86400)
  ## Look for sustained elevated bands.
  if (length(panel2$rate) >= 5L) {
    elev <- panel2$rate > 0.05    # 5% in the window
    n_elev_runs <- length(rle(elev)$lengths[rle(elev)$values])
    if (any(elev) && max(rle(elev)$lengths[rle(elev)$values], na.rm = TRUE) >= 3L) {
      notes <- c(notes,
        "Panel 2: sustained band of elevated flag rate detected -- this is the migration-over-flagging signature. Consider filtering that window or running it through state-conditional analysis.")
    }
  }

  ## Panel 3: per-class breakdown (class-aware taxonomy)
  ## Counts the canonical error_class strings produced by mt_clean_track().
  ec <- xf$error_class
  ec <- ec[!is.na(ec)]
  ## Enumerate the full vocabulary so the bar plot has a stable order.
  vocab <- c("consensus", "geometric_spike", "state_anomaly",
             "kinematic_confluence", "block", "physiological",
             "state_transition_buffered")
  combos <- vapply(vocab, function(v) sum(ec == v), integer(1))
  combos <- combos[combos > 0]
  if (sum(combos) > 0) {
    n_high_conf <- sum(combos[intersect(names(combos),
                          c("consensus", "geometric_spike",
                            "state_anomaly", "kinematic_confluence",
                            "block", "physiological"))])
    high_conf_frac <- n_high_conf / sum(combos)
    if (high_conf_frac < 0.5) {
      dominant <- names(combos)[which.max(combos)]
      notes <- c(notes, sprintf(
        "Panel 3: only %.0f%% of flags fall in high-confidence multi-rule classes; dominant class is '%s'. Inspect the per-detector columns directly to see what is driving the flags.",
        100 * high_conf_frac, dominant))
    }
  }

  ## Panel 4: cumulative by iteration
  iters <- xf$flag_iteration
  iters <- iters[!is.na(iters)]
  max_iter <- if (length(iters)) max(iters, na.rm = TRUE) else 0L
  cum_per_iter <- if (length(iters)) {
    vapply(seq_len(max_iter), function(k) sum(iters <= k), integer(1))
  } else integer(0)
  conv <- attr(x, "convergence")
  if (!is.null(conv) && conv == "flag_fraction_exceeded") {
    notes <- c(notes,
      "Panel 4: iteration hit max_flag_fraction abort -- the detectors are not converging on a stable flag set. See Panels 1 and 2 for state-bimodality / migration signatures.")
  }

  ## Panel 5: run lengths
  run_lengths <- if (any(is_out)) {
    rle_out <- rle(as.integer(is_out))
    rle_out$lengths[rle_out$values == 1L]
  } else integer(0)
  if (length(run_lengths) > 0L) {
    long_runs <- sum(run_lengths >= 10L)
    if (long_runs > 0L &&
        long_runs * 10L > 0.5 * sum(run_lengths)) {
      notes <- c(notes, sprintf(
        "Panel 5: %d runs of >=10 consecutive flagged fixes -- a sustained block-shaped flag pattern. Either supply v_max so block expansion can resolve genuine error clusters, or check Panel 2 for behavioural-state mismatch.",
        long_runs))
    }
  }

  ## Panel 6: per-individual rates -- only if multi-track
  high_rate_inds <- by_individual[by_individual$pct > 2, ]
  if (is_multi && nrow(high_rate_inds) > 0L) {
    notes <- c(notes, sprintf(
      "Panel 6: %d / %d individual(s) flagged at >2%% (likely state-bimodal): %s. See state-conditional vignette (in preparation).",
      nrow(high_rate_inds), nrow(by_individual),
      paste(substr(high_rate_inds$individual, 1, 30),
            collapse = ", ")))
  }

  ## --- render -------------------------------------------------------
  if (plot) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    n_panels <- if (is_multi) 6L else 5L
    nr <- if (n_panels == 6L) 2L else 2L
    nc <- if (n_panels == 6L) 3L else 3L     # 2x3 layout fits 5 or 6 panels
    ## The mfrow layout auto-shrinks base text; counteract it so the
    ## titles, labels and tick marks stay legible when saved to a file
    ## (Elisa feedback 2026-06-01).  cex_scale lets the user push it
    ## further.
    graphics::par(mfrow = c(nr, nc),
                  mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0),
                  cex.main = cex_scale, cex.lab = cex_scale,
                  cex.axis = 0.9 * cex_scale)

    .panel_speed_density(v_pos, v_max_focus, modes, cex_scale)
    .panel_flag_rate_time(panel2, cex_scale)
    .panel_detector_contribution(combos, cex_scale)
    .panel_iteration_cumulative(cum_per_iter, conv,
                                  attr(x, "v_max_used"), cex_scale)
    .panel_run_lengths(run_lengths, cex_scale)
    if (is_multi) .panel_per_individual(by_individual, cex_scale)

    graphics::mtext(sprintf(
      "mt_clean_track diagnostics -- focus: %s (%d / %d flagged, %.3f%%)",
      individual, sum(is_out), length(is_out),
      100 * mean(is_out)),
      outer = TRUE, line = 1, cex = 1.05 * cex_scale, font = 2)
  }

  ## --- print interpretive notes ------------------------------------
  if (length(notes) > 0L) {
    say("=== mt_diagnose_clean_track: concerns flagged ===")
    for (n in notes) say(n)
  } else {
    say("=== mt_diagnose_clean_track: no concerns flagged. ===")
  }

  invisible(list(
    individual    = individual,
    by_individual = by_individual,
    run_lengths   = run_lengths,
    modes         = modes,
    notes         = notes
  ))
}


# ---- helpers --------------------------------------------------------------

## Find substantive modes in log(speed) using the same definition as
## the auto-cap mode-position gate: basin holds >= 2% of fixes AND
## peak density >= 5% of global peak.
##
## @keywords internal
.find_substantive_modes <- function(v_pos,
                                      mode_min_frac    = 0.02,
                                      mode_min_density = 0.05,
                                      n_grid = 512) {
  if (length(v_pos) < 30L) return(numeric(0))
  log_s <- log(v_pos)
  d <- stats::density(log_s, n = n_grid)
  dy <- diff(d$y)
  is_max <- which(c(FALSE, dy[-length(dy)] > 0 & dy[-1L] < 0, FALSE))
  is_min <- which(c(FALSE, dy[-length(dy)] < 0 & dy[-1L] > 0, FALSE))
  if (length(is_max) == 0L) return(numeric(0))
  global_peak <- max(d$y[is_max])
  modes <- numeric(0)
  for (idx in is_max) {
    left  <- is_min[is_min < idx]
    right <- is_min[is_min > idx]
    lo <- if (length(left))  d$x[max(left)]   else -Inf
    hi <- if (length(right)) d$x[min(right)]  else  Inf
    frac <- mean(log_s >= lo & log_s <= hi)
    dr   <- d$y[idx] / global_peak
    if (frac >= mode_min_frac && dr >= mode_min_density) {
      modes <- c(modes, exp(d$x[idx]))
    }
  }
  modes
}


## Rolling flag rate computed at the midpoints of consecutive
## non-overlapping windows of size `window_secs`.
##
## @keywords internal
.rolling_flag_rate <- function(ts, is_out, window_secs, min_windows = 8L) {
  if (length(ts) < 2L) {
    return(list(time = numeric(0), rate = numeric(0),
                window_secs = window_secs))
  }
  t_min <- min(ts, na.rm = TRUE)
  t_max <- max(ts, na.rm = TRUE)
  if (!is.finite(t_min) || !is.finite(t_max)) {
    return(list(time = numeric(0), rate = numeric(0),
                window_secs = window_secs))
  }
  ## Auto-shrink for short tracks: the requested window (default
  ## 7 days) yields < min_windows panels on tracks spanning only a few
  ## days, which previously fell through to a "too few windows"
  ## placeholder.  Target at least min_windows windows by narrowing the
  ## width to span / min_windows when the track is too short for the
  ## requested width.  Long tracks keep the requested width unchanged.
  span_s <- t_max - t_min
  if (is.finite(span_s) && span_s > 0 &&
      span_s / window_secs < min_windows) {
    window_secs <- span_s / min_windows
  }
  breaks <- seq(t_min, t_max, by = window_secs)
  if (length(breaks) < 2L) breaks <- c(t_min, t_max)
  midpts <- (utils::head(breaks, -1L) + utils::tail(breaks, -1L)) / 2
  rate   <- vapply(seq_along(midpts), function(k) {
    sel <- ts >= breaks[k] & ts < breaks[k + 1L]
    if (!any(sel)) NA_real_ else mean(is_out[sel], na.rm = TRUE)
  }, numeric(1))
  list(time = midpts, rate = rate, window_secs = window_secs)
}


## --- panel renderers ------------------------------------------------------

.panel_speed_density <- function(v_pos, v_max, modes, cex_scale = 1) {
  if (length(v_pos) < 10L) {
    graphics::plot.new(); graphics::title("1. step-speed density (log scale)")
    graphics::text(0.5, 0.5, "too few positive speeds",
                   cex = 1.05 * cex_scale)
    return(invisible(NULL))
  }
  log_s <- log(v_pos)
  d <- stats::density(log_s, n = 512)
  graphics::plot(d$x, d$y, type = "l", col = "grey30", lwd = 1.4,
                 xlab = "step speed (m/s, log scale)", ylab = "density",
                 main = "1. step-speed density (log scale)")
  ## annotate substantive modes
  if (length(modes) > 0L) {
    graphics::abline(v = log(modes), col = "steelblue4", lty = 2, lwd = 1.4)
    graphics::text(log(modes), max(d$y) * 0.95,
                   sprintf("%.2f m/s", modes),
                   col = "steelblue4", cex = 0.85 * cex_scale,
                   pos = 4, srt = 0)
  }
  if (is.finite(v_max)) {
    graphics::abline(v = log(v_max), col = "firebrick", lty = 1, lwd = 2)
    graphics::text(log(v_max), max(d$y) * 0.5,
                   sprintf("v_max = %.1f m/s", v_max),
                   col = "firebrick", cex = 0.9 * cex_scale, pos = 4, srt = 90)
  }
  graphics::legend("topleft",
                   legend = c("density",
                              sprintf("substantive mode (n=%d)", length(modes)),
                              "v_max applied"),
                   col = c("grey30", "steelblue4", "firebrick"),
                   lty = c(1, 2, 1), lwd = c(1.4, 1.4, 2),
                   bty = "n", cex = 0.85 * cex_scale)
}

.panel_flag_rate_time <- function(panel2, cex_scale = 1) {
  if (length(panel2$rate) < 2L) {
    graphics::plot.new(); graphics::title("2. flag rate vs time")
    graphics::text(0.5, 0.55, "too few windows", cex = 1.1 * cex_scale)
    graphics::text(0.5, 0.38,
                   "track spans too little time to resolve a flag-rate series",
                   cex = 0.85 * cex_scale, col = "grey40")
    return(invisible(NULL))
  }
  ## convert seconds to POSIXct for axis
  tt <- as.POSIXct(panel2$time, origin = "1970-01-01", tz = "UTC")
  graphics::plot(tt, 100 * panel2$rate, type = "l", lwd = 1.4,
                 col = "grey30",
                 xlab = "time", ylab = "flag rate (% in window)",
                 main = "2. flag rate vs time")
  graphics::abline(h = c(0.5, 2, 5), col = c("grey60", "orange", "firebrick"),
                   lty = 2, lwd = 1)
  ## annotate the window width actually used (may be auto-shrunk on
  ## short tracks)
  win_secs <- panel2$window_secs
  if (!is.null(win_secs) && is.finite(win_secs)) {
    win_lab <- if (win_secs >= 86400)
      sprintf("window = %.1f d", win_secs / 86400)
    else
      sprintf("window = %.1f h", win_secs / 3600)
    graphics::mtext(win_lab, side = 3, line = 0.2,
                    cex = 0.8 * cex_scale, col = "grey40")
  }
  graphics::legend("topright",
                   title = "rate thresholds",
                   legend = c("0.5% (healthy)", "2% (concerning)",
                              "5% (alarming)"),
                   col = c("grey60", "orange", "firebrick"),
                   lty = 2, bty = "n", cex = 0.85 * cex_scale)
}

.panel_detector_contribution <- function(combos, cex_scale = 1) {
  if (sum(combos) == 0L) {
    graphics::plot.new(); graphics::title("3. per-detector activity")
    graphics::text(0.5, 0.5, "no flags", cex = 1.05 * cex_scale)
    return(invisible(NULL))
  }
  cols <- c(
    "consensus"                  = "forestgreen",
    "geometric_spike"            = "purple",
    "state_anomaly"              = "mediumvioletred",
    "kinematic_confluence"       = "tan4",
    "block"                      = "darkred",
    "physiological"              = "grey40",
    "state_transition_buffered"  = "yellow")
  ord <- order(combos, decreasing = FALSE)   # smallest at bottom for horiz
  ## widen left margin so the long class names fit
  op2 <- graphics::par(mar = c(4, 9, 3, 2))
  on.exit(graphics::par(op2), add = TRUE)
  bp <- graphics::barplot(combos[ord],
                          col = cols[names(combos)[ord]],
                          horiz = TRUE, las = 1,
                          main = "3. error class breakdown",
                          xlab = "fixes (cumulative across iterations)",
                          cex.names = 0.85 * cex_scale,
                          xlim = c(0, max(combos) * 1.15))
  graphics::text(combos[ord], bp, labels = combos[ord],
                 pos = 4, cex = 0.85 * cex_scale, xpd = NA)
}

.panel_iteration_cumulative <- function(cum_per_iter, conv, v_max_used,
                                          cex_scale = 1) {
  if (length(cum_per_iter) == 0L) {
    graphics::plot.new(); graphics::title("4. cumulative flags by iteration")
    graphics::text(0.5, 0.5, "no flags -- nothing to plot",
                   cex = 1.05 * cex_scale)
    return(invisible(NULL))
  }
  iters <- seq_along(cum_per_iter)
  if (length(iters) == 1L) {
    ## Single iteration: the cascade saturated on the first pass.  A
    ## lone full-width bar reads as an uninformative solid box, so draw
    ## a compact text summary instead (Kamran/Elisa feedback 2026-06-01).
    n_flag <- cum_per_iter[1L]
    graphics::plot.new(); graphics::title("4. cumulative flags by iteration")
    graphics::text(0.5, 0.62,
                   sprintf("%d fix%s flagged on the first pass",
                           n_flag, if (n_flag == 1L) "" else "es"),
                   cex = 1.1 * cex_scale)
    graphics::text(0.5, 0.40, "converged in 1 iteration",
                   cex = 1.0 * cex_scale, col = "forestgreen", font = 2)
    graphics::text(0.5, 0.25,
                   "(no second pass needed -- nothing to accumulate)",
                   cex = 0.85 * cex_scale, col = "grey40")
    return(invisible(NULL))
  }
  graphics::plot(iters, cum_per_iter, type = "b", pch = 19, lwd = 1.4,
                 col = "grey20",
                 xlab = "iteration", ylab = "cumulative flags",
                 main = "4. cumulative flags by iteration")
  if (!is.null(conv)) {
    sub <- sprintf("convergence: %s", conv)
    if (conv == "flag_fraction_exceeded") sub <- paste(sub, "(ABORT)")
    graphics::mtext(sub, side = 3, line = 0.2, cex = 0.85 * cex_scale,
                    col = if (conv == "flag_fraction_exceeded") "firebrick"
                          else "grey30")
  }
}

.panel_run_lengths <- function(run_lengths, cex_scale = 1) {
  if (length(run_lengths) == 0L) {
    graphics::plot.new(); graphics::title("5. consecutive-flag run lengths")
    graphics::text(0.5, 0.5, "no flags", cex = 1.05 * cex_scale)
    return(invisible(NULL))
  }
  ## bin: 1, 2, 3, 4, 5-9, 10+
  bins <- c(1, 2, 3, 4, 5, 10, Inf)
  labels <- c("1", "2", "3", "4", "5-9", "10+")
  bin_idx <- findInterval(run_lengths, bins, all.inside = TRUE)
  counts <- tabulate(bin_idx, nbins = length(labels))
  cols <- c("grey80", "grey65", "grey50", "grey40", "orange", "firebrick")
  bp <- graphics::barplot(counts, names.arg = labels, col = cols,
                          main = "5. consecutive-flag run lengths",
                          xlab = "run length (consecutive flagged fixes)",
                          ylab = "count of runs",
                          cex.names = 0.95 * cex_scale)
  graphics::text(bp, counts, labels = counts, pos = 3,
                 cex = 0.9 * cex_scale, xpd = NA)
}

.panel_per_individual <- function(by_individual, cex_scale = 1) {
  if (nrow(by_individual) == 0L) {
    graphics::plot.new(); graphics::title("6. per-individual flag rates")
    return(invisible(NULL))
  }
  d <- by_individual
  d <- d[order(-d$pct), ]
  n <- nrow(d)
  cols <- ifelse(d$pct > 5, "firebrick",
            ifelse(d$pct > 2, "orange",
              ifelse(d$pct > 0.5, "tan", "forestgreen")))
  ## widen left margin a touch so individual names fit
  op2 <- graphics::par(mar = c(4, 8, 3, 1))
  on.exit(graphics::par(op2), add = TRUE)
  graphics::plot(d$pct, seq_len(n), pch = 19, col = cols, cex = 1.2 * cex_scale,
                 yaxt = "n", xlim = c(0, max(d$pct, 5) * 1.05),
                 xlab = "flagged (%)", ylab = "",
                 main = "6. per-individual flag rates")
  graphics::axis(2, at = seq_len(n),
                 labels = substr(d$individual, 1, 22),
                 las = 1, cex.axis = 0.85 * cex_scale)
  graphics::abline(v = c(0.5, 2, 5),
                   col = c("grey60", "orange", "firebrick"),
                   lty = 2)
}

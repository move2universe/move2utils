#' Sequential outlier detection for movement data
#'
#' Detects outliers by walking through a track from confirmed-good
#' locations, evaluating each candidate step against the same joint
#' probability (Equation 1) used by \code{\link{mt_flag_outliers}}.
#'
#' Three scanning strategies are available:
#' \describe{
#'   \item{\code{"forward-backward"}}{Scan forward from the first
#'     location and backward from the last; combine via max probability.}
#'   \item{\code{"greedy"}}{Start from the location with the highest
#'     simultaneous joint probability and expand outward in both
#'     directions.  Requires running \code{mt_flag_outliers} first.}
#'   \item{\code{"random"}}{Launch scans from multiple random anchor
#'     points; each location's probability is the median across scans.}
#' }
#'
#' @details
#' \strong{Relationship to the four-primitive cascade.}  The unified
#' cleaner \code{\link{mt_clean_track}} uses \code{simultaneous}
#' thresholding on the joint-probability surface (via
#' \code{\link{mt_flag_outliers}}) and combines it with three other
#' detectors (bridge, detour, speed-cap) under a class-aware flag
#' rule.  \code{mt_sequential_outliers} is an alternative
#' \emph{strategy} on the same probability surface: rather than
#' applying a single threshold to all fixes at once, it walks from
#' confirmed-good locations and evaluates each step as a transition.
#' This catches block-shaped errors whose interior steps look
#' ordinary in transitions (the simultaneous threshold misses them
#' because each interior step is locally plausible).  The cascade
#' addresses block-shaped contamination via bridge + detour + the
#' topological block-expansion step, so \code{mt_clean_track} is the
#' recommended default.  Reach for \code{mt_sequential_outliers}
#' when (i) you want an orthogonal cross-check on the cascade
#' results, (ii) you have a tightly-clustered burst of contamination
#' on an otherwise short track where the cascade's iteration loop
#' converges slowly, or (iii) you are calibrating thresholds on the
#' joint-probability surface specifically and want to compare
#' simultaneous vs sequential evaluation.
#'
#' @param x A \code{move2} object.
#' @param reference Optional \code{move2} object from which the
#'   probability surfaces are built.  If NULL, built from \code{x}.
#' @param scan Scanning strategy: \code{"forward-backward"} (default),
#'   \code{"greedy"}, or \code{"random"}.
#' @param n_random Number of random anchor points for
#'   \code{scan = "random"}.  Default: 10.
#' @param autodiff_alpha Exponent on the auto-difference terms
#'   (Equation 1).  Default: 0.5.
#' @param threshold Minimum joint probability for accepting a step.
#'   If NULL (default), set to the 0.5th percentile of the reference
#'   joint probability distribution.
#' @param max_skip Maximum consecutive locations to skip before
#'   force-advancing.  Default: 100.
#' @param time_normalize Logical; if TRUE (default), use speed and
#'   angular velocity.
#' @param plot Logical; if TRUE, plot the results.
#' @param anchor_corruption_threshold Numeric in \code{(0, 1)}.  If a
#'   track's per-track flag rate exceeds this value the function
#'   emits a warning naming the anchor-corruption signature.  Default
#'   \code{0.30}: above 30\% of fixes flagged from a single-anchor
#'   scan, the anchor is the more likely failure mode than the data.
#'   Set to \code{NULL} to disable.
#'
#' @section Anchor-corruption diagnostic:
#' Sequential scanning evaluates each step relative to a confirmed-
#' good anchor.  When the anchor is itself corrupt (a tag-deployment
#' glitch on the first or last fix; an Argos bias spike that
#' coincides with the strategy's chosen anchor in
#' \code{scan = "greedy"}; or all anchors in \code{scan = "random"}
#' falling within a contaminated cluster), the scan immediately
#' classifies every legitimate fix as anomalous because every
#' legitimate step looks "anomalous" relative to the corrupt anchor.
#' The \code{max_skip} safety net catches the runaway case eventually,
#' but below that threshold the failure mode is silent.
#'
#' This diagnostic surfaces the failure mode by warning when more
#' than \code{anchor_corruption_threshold} of a track's fixes are
#' flagged.  Real outlier rates are typically <5\%; rates above 30\%
#' indicate the scan is reflecting the anchor's perspective rather
#' than the data's.  Recovery: try \code{scan = "forward-backward"}
#' (already the default; combines both endpoints), supply a clean
#' \code{reference =} object, or use \code{\link{mt_clean_track}}
#' which is anchor-free.
#'
#' @section Two-threshold calibration (fixed 2026-05-12):
#' The scan-time scoring formula has two regimes: a partial-formula
#' first step from any anchor (\code{prob = stp}, no autodifference
#' contribution because there is no prior anchor or post-\code{max_skip}
#' reset wiped \code{prev_v}/\code{prev_w}) and a full-formula in-scan
#' step (\code{prob = stp * (dsp * dtp)^autodiff_alpha}).  The
#' algorithm now calibrates two thresholds (\code{threshold_partial}
#' and \code{threshold_full}) -- each the 0.5\%-ile of its matching
#' reference distribution -- and the scan compares each step's score
#' to its regime-matched threshold.  Pre-fix (before 2026-05-12) the
#' algorithm used a single threshold calibrated on the full formula;
#' on multi-state tracks where the autodifference KDE returned very
#' large densities at the bimodal Delta-step peaks, the full-formula
#' threshold sat orders of magnitude above typical partial-formula
#' first-step scores and the scan flagged ~all fixes (WH17 over-flag
#' at 99.99\%).  Post-fix on WH17 the flag rate drops to ~1.4\%; on
#' single-state synthetic CPF tracks the two thresholds happen to
#' agree numerically and the per-fix flag set is byte-identical.
#' User-supplied \code{threshold = X} sets both thresholds to
#' \code{X} (back-compat).
#'
#' The 0.5\%-percentile heuristic remains the calibration target;
#' a per-fix p-value-based threshold (parametric null on the
#' reference log-probability) is the natural improvement direction
#' but is deferred (see \code{FUTURE_IMPROVEMENTS.md} M-12).
#'
#' Documented workflow recommendations remain unchanged: on heavily
#' contaminated multi-state tracks, supply a clean
#' \code{reference =} object, or use \code{\link{mt_flag_outliers}}
#' for initial screening, or use \code{\link{mt_clean_track}}
#' (anchor-free).  Deriving a state column from speed and
#' dispatching per-state does NOT help -- it fragments the failure
#' AND confounds outlier signal with state signal on synthetic
#' block-shape contamination (CPF_D).
#'
#' @return The input \code{move2} object with added columns:
#'   \code{is_outlier} and \code{seq_joint_prob}.
#'
#' @references
#' Safi, K. (in preparation). Self-thresholding hierarchical
#' outlier-detection for animal movement tracks. Companion paper to
#' the \pkg{move2utils} R package. Preprint: bioRxiv (DOI forthcoming).
#'
#' @seealso \code{\link{mt_clean_track}} (recommended unified
#'   cleaner); \code{\link{mt_flag_outliers}} (probability primitive
#'   with simultaneous thresholding -- the basis this function
#'   scans over); \code{\link{mt_combined_outliers}}
#'   (majority vote across simultaneous + sequential strategies);
#'   \code{\link{mt_persistence_score}} (multi-scale persistence
#'   annotation -- works on the output of any flagger including
#'   this one).
#'
#' @examples
#' \dontrun{
#' ## Forward-backward scan from the track endpoints:
#' res <- mt_sequential_outliers(track, scan = "forward-backward")
#' summary(res$seq_joint_prob)
#' }
#'
#' @export
mt_sequential_outliers <- function(x, reference = NULL,
                                   scan = "forward-backward",
                                   n_random = 10,
                                   autodiff_alpha = 0.5,
                                   threshold = NULL,
                                   max_skip = 100,
                                   time_normalize = TRUE,
                                   anchor_corruption_threshold = 0.30,
                                   plot = FALSE) {
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  .reject_empty_geometry(x, "mt_sequential_outliers()")
  scan <- match.arg(scan, c("forward-backward", "greedy", "random"))
  if (!is.null(anchor_corruption_threshold)) {
    if (!is.numeric(anchor_corruption_threshold) ||
        length(anchor_corruption_threshold) != 1L ||
        is.na(anchor_corruption_threshold) ||
        anchor_corruption_threshold <= 0 ||
        anchor_corruption_threshold >= 1) {
      rlang::abort("`anchor_corruption_threshold` must be NULL or a scalar in (0, 1).",
                   class = "move2utils_mt_sequential_outliers_bad_anchor_corruption_threshold")
    }
  }

  ## Stratified projection handling: project once to a local metric CRS so
  ## every downstream geometric computation is in metres.  A move2 object
  ## carries its CRS, so this is done internally rather than forcing the
  ## caller to project first.
  orig_crs <- sf::st_crs(x)
  x <- .to_canonical_aeqd(x)
  if (!is.null(reference)) reference <- .to_canonical_aeqd(reference)

  ## ---- build probability surfaces ----
  ref <- if (!is.null(reference)) reference else x
  message("Building probability surfaces...")

  ref_step <- as.numeric(move2::mt_distance(ref, units = "m"))
  ## projection-agnostic, track-aware turn angles (Cartesian on projected
  ## input, geographic on longlat); move2::mt_turnangle() is longlat-only
  ## and errors on a metric CRS, which broke this function on projected
  ## tracks.  Step lengths via mt_distance() are already metric-correct
  ## (metres) in either CRS.
  ref_turn <- as.numeric(.turn_angles_fast(ref))
  if (time_normalize) {
    ref_tl <- as.numeric(move2::mt_time_lags(ref))
    ref_na <- is.na(ref_tl)
    ref_step <- ifelse(ref_na, NA_real_, ref_step / ref_tl)
    ref_turn <- ifelse(ref_na, NA_real_, ref_turn / ref_tl)
  }
  ref_deltaStep <- c(NA, diff(ref_step))
  ref_deltaTurn <- c(NA, diff(ref_turn))

  ## gap lengths for autodifference scaling
  if (time_normalize) {
    ref_tl_num <- as.numeric(move2::mt_time_lags(ref))
    ref_delta_gaps <- c(NA, (ref_tl_num[-length(ref_tl_num)] +
                             ref_tl_num[-1]) / 2)
  } else {
    ref_delta_gaps <- rep(1, length(ref_step))
  }

  hist2d <- .turn_step_hist(ref_turn, ref_step)

  ## pre-extract the 2D histogram into a matrix for fast lookup
  ## (avoids per-step terra::extract S4 overhead in the scan loop)
  hist_mat <- as.matrix(hist2d, wide = TRUE)
  hist_ext <- as.vector(terra::ext(hist2d))  # xmin, xmax, ymin, ymax
  hist_res <- terra::res(hist2d)              # dx, dy
  hist_nr  <- nrow(hist2d)
  hist_nc  <- ncol(hist2d)

  ## fast 2D histogram lookup: (angular_velocity, speed) → density
  hist_lookup <- function(w, v) {
    col <- as.integer(floor((w - hist_ext[1]) / hist_res[1]) + 1)
    row <- as.integer(hist_nr - floor((v - hist_ext[3]) / hist_res[2]))
    col <- pmin(pmax(col, 1L), hist_nc)
    row <- pmin(pmax(row, 1L), hist_nr)
    hist_mat[cbind(row, col)]
  }
  ga_step <- .gap_aware_autodiff(ref_deltaStep, ref_delta_gaps)
  ga_turn <- .gap_aware_autodiff(ref_deltaTurn, ref_delta_gaps)

  ## ---- determine thresholds (two-threshold fix, 2026-05-12) ---------
  ##
  ## The scan-time scoring formula has two regimes:
  ##   (a) FIRST step from any anchor (prev_anchor is NA): the
  ##       autodifference is undefined, so score_triplet returns
  ##         prob = stp                  (the 2D-histogram-only score)
  ##   (b) IN-scan steps (prev_anchor is set): the full formula
  ##         prob = stp * (dsp * dtp)^autodiff_alpha
  ##
  ## Pre-fix: a single threshold was calibrated as the 0.5%-ile of
  ## the FULL formula on the reference data and applied to BOTH
  ## regimes.  On multi-state tracks (e.g.\ a flying/perching gull
  ## like WH17) the gap-aware KDE on the bimodal Delta-step
  ## distribution produces very large densities; the full-formula
  ## threshold sits orders of magnitude above the partial-formula
  ## value of typical first-step scores, so every first-step is
  ## flagged.  The scan's anchor never advances past its initial
  ## point, max_skip triggers, and ~all fixes get flagged ("anchor-
  ## corruption" signature -- but driven by scale mismatch, not a
  ## corrupt anchor).  Empirically on WH17: 99.99% flagged.  See
  ## FUTURE_IMPROVEMENTS.md "Item E" for the full diagnosis.
  ##
  ## Fix: calibrate TWO thresholds, one per regime, each on its
  ## matching reference distribution.  Each scan-time score is then
  ## compared to a same-formula threshold and the scale mismatch is
  ## eliminated by construction.  Single-state tracks where the two
  ## thresholds happen to agree (e.g.\ CPF synthetic where
  ## (dsp*dtp)^alpha ~ 1 because the autodiff KDE returns ~1 at
  ## the bulk of a unimodal Delta-step distribution) are byte-
  ## identical to the pre-fix behaviour.
  ##
  ## User-supplied `threshold = X` (a single scalar) sets BOTH
  ## thresholds to X, preserving the pre-fix override behaviour.
  ref_stp <- terra::extract(hist2d, cbind(ref_turn, ref_step))[, 1]
  ref_stp <- pmax(ref_stp, .Machine$double.xmin, na.rm = TRUE)
  n_ref <- length(ref_step)
  ref_dsp <- rep(1, n_ref); ref_dtp <- rep(1, n_ref)
  ok <- !is.na(ref_deltaStep) & !is.na(ref_deltaTurn) &
        !is.na(ref_delta_gaps) & ref_delta_gaps > 0
  if (any(ok)) {
    ## vectorised: one scale / KDE call per branch, not per row
    gaps_ok <- ref_delta_gaps[ok]
    s_s <- ga_step$scale_fun(gaps_ok)
    pos_s <- !is.na(s_s) & s_s > 0
    if (any(pos_s)) {
      idx <- which(ok)[pos_s]
      s_vec <- s_s[pos_s]
      ref_dsp[idx] <- pmax(
        ga_step$kde_fun(ref_deltaStep[idx] / s_vec) / s_vec,
        .Machine$double.xmin)
    }
    s_t <- ga_turn$scale_fun(gaps_ok)
    pos_t <- !is.na(s_t) & s_t > 0
    if (any(pos_t)) {
      idx <- which(ok)[pos_t]
      s_vec <- s_t[pos_t]
      ref_dtp[idx] <- pmax(
        ga_turn$kde_fun(.wrap_angle(ref_deltaTurn[idx]) / s_vec) / s_vec,
        .Machine$double.xmin)
    }
  }
  ref_joint <- ref_stp * (ref_dsp * ref_dtp)^autodiff_alpha
  if (is.null(threshold)) {
    ## Calibrate the thresholds from the REAL densities only, excluding the
    ## .Machine$double.xmin floor that scores collapse to when a fix's
    ## (turn, step) falls outside the 2-D histogram support.  Those floored
    ## fixes are the maximally-anomalous ones (the outliers themselves); if
    ## more than 0.5% of the reference is floored -- which happens on short
    ## tracks with a few large spikes -- the raw 0.5%ile lands on the floor,
    ## the threshold degenerates to double.xmin, and "score < threshold"
    ## can never fire, so nothing is flagged.  Excluding the floor makes the
    ## threshold a real low-density value that the floored outliers fall
    ## below.  (Latent bug surfaced when projection canonicalisation made
    ## CPF_C's histogram sparse; see mt_sequential_outliers tests.)
    floor_v       <- .Machine$double.xmin
    ref_stp_pos   <- ref_stp[!is.na(ref_stp)     & ref_stp     > floor_v]
    ref_joint_pos <- ref_joint[!is.na(ref_joint) & ref_joint   > floor_v]
    threshold_partial <- if (length(ref_stp_pos))
      stats::quantile(ref_stp_pos,   0.005, na.rm = TRUE) else floor_v
    threshold_full    <- if (length(ref_joint_pos))
      stats::quantile(ref_joint_pos, 0.005, na.rm = TRUE) else floor_v
    ## Back-compat alias: a single scalar consumers may inspect.  The
    ## actual scan uses the regime-matched threshold below.  Use the
    ## full-formula one (the pre-fix definition) so attr-consumers
    ## are not surprised.
    threshold <- threshold_full
  } else {
    ## User-supplied scalar: apply to both regimes (back-compat).
    threshold_partial <- threshold
    threshold_full    <- threshold
  }
  message(sprintf(
    "Sequential scan (%s): threshold_partial = %.4g, threshold_full = %.4g",
    scan, threshold_partial, threshold_full))

  ## ---- determine time unit ----
  tl_raw <- move2::mt_time_lags(ref)
  time_unit <- if (inherits(tl_raw, "units")) units(tl_raw) else NULL

  ## ---- greedy: get simultaneous probs for anchor selection ----
  sim_probs <- NULL
  if (scan == "greedy") {
    message("Running simultaneous method to find best anchor...")
    r_sim <- suppressMessages(
      mt_flag_outliers(x, plot = FALSE, time_normalize = time_normalize)
    )
    sim_probs <- r_sim$joint_prob
  }

  ## ---- per-individual processing ----
  ids <- move2::mt_track_id(x)
  unique_ids <- unique(ids)
  all_outlier <- rep(FALSE, nrow(x))
  all_prob <- rep(NA_real_, nrow(x))

  for (uid in unique_ids) {
    idx <- which(ids == uid)
    xi <- x[idx, ]
    ni <- length(idx)
    if (ni < 3) next

    ## pre-compute metric coordinates using LOCAL latitude for the
    ## Coordinates are metres: canonicalisation to a local metric AEQD is
    ## handled once at the entry point (see .to_canonical_aeqd() at the top),
    ## so every downstream computation here simply uses them directly.
    coords <- sf::st_coordinates(xi)
    cx <- coords[, 1]
    cy <- coords[, 2]

    posix_t <- move2::mt_time(xi)
    if (!is.null(time_unit)) {
      dt_secs <- as.numeric(difftime(posix_t, posix_t[1], units = "secs"))
      times <- as.numeric(
        units::set_units(units::set_units(dt_secs, "s"),
                          time_unit, mode = "standard"))
    } else {
      times <- as.numeric(posix_t)
    }

    ## ---- score_triplet: Equation 1 for a (prev, anchor, cand) step ----
    ## prev_dt is the time lag of the previous accepted step (for gap-
    ## aware autodifference scaling).
    ##
    ## Returns is_partial = TRUE when the autodifference contribution
    ## is undefined (no prior anchor, or no prior step speed/turn after
    ## a max_skip reset) -- in that case prob = stp only, scoped against
    ## threshold_partial.  Otherwise prob is the full formula and
    ## scoped against threshold_full.  See the two-threshold rationale
    ## above for why this distinction is necessary.
    score_triplet <- function(prev, anchor, cand, prev_v, prev_w, prev_dt) {
      dx <- cx[cand] - cx[anchor]; dy <- cy[cand] - cy[anchor]
      dist_m <- sqrt(dx * dx + dy * dy)
      dt <- abs(times[cand] - times[anchor])
      if (is.na(dt) || dt <= 0)
        return(list(prob = NA, v = NA, w = NA, dt = NA, is_partial = NA))

      v <- if (time_normalize) dist_m / dt else dist_m

      if (is.na(prev)) {
        w <- 0
        stp <- max(hist_lookup(w, v), .Machine$double.xmin, na.rm = TRUE)
        return(list(prob = stp, v = v, w = NA_real_, dt = dt,
                    is_partial = TRUE))
      }

      dx1 <- cx[anchor] - cx[prev];  dy1 <- cy[anchor] - cy[prev]
      ## dx, dy already computed above for distance
      angle <- atan2(dx * dy1 - dy * dx1, dx * dx1 + dy * dy1)
      w <- if (time_normalize) angle / dt else angle

      stp <- max(hist_lookup(w, v), .Machine$double.xmin, na.rm = TRUE)

      ## gap-aware autodifferences
      dv <- if (!is.na(prev_v)) v - prev_v else NA_real_
      dw <- if (!is.na(prev_w)) w - prev_w else NA_real_

      ## the gap for the autodifference is the average of the current
      ## and previous step durations (time between step midpoints)
      delta_gap <- if (!is.na(prev_dt) && prev_dt > 0) {
        (prev_dt + dt) / 2
      } else {
        dt
      }

      dsp <- 1; dtp <- 1
      has_dsp <- FALSE; has_dtp <- FALSE
      if (!is.na(dv)) {
        s <- ga_step$scale_fun(delta_gap)
        if (!is.na(s) && s > 0) {
          dsp <- max(ga_step$kde_fun(dv / s) / s, .Machine$double.xmin)
          has_dsp <- TRUE
        }
      }
      if (!is.na(dw)) {
        s <- ga_turn$scale_fun(delta_gap)
        if (!is.na(s) && s > 0) {
          dtp <- max(ga_turn$kde_fun(.wrap_angle(dw) / s) / s,
                      .Machine$double.xmin)
          has_dtp <- TRUE
        }
      }
      ## Both default 1 (no autodiff contribution) -> equivalent to
      ## the partial-formula case numerically.
      is_partial <- !has_dsp && !has_dtp

      list(prob = stp * (dsp * dtp)^autodiff_alpha, v = v, w = w, dt = dt,
           is_partial = is_partial)
    }

    ## ---- single-direction scan from a starting anchor ----
    ## direction = +1 (forward) or -1 (backward)
    ##
    ## Returns a list(prob, is_partial) -- both length-ni vectors -- so
    ## the wrapper can pick the regime-matched threshold per fix.  The
    ## scan's accept-or-skip decision uses the SAME regime-matched
    ## threshold (Formula A -> threshold_partial, Formula B -> threshold_full).
    run_scan <- function(start, direction) {
      probs   <- rep(NA_real_, ni)
      partial <- rep(NA, ni)
      anchor <- start
      prev_anchor <- NA_integer_
      prev_v <- NA_real_; prev_w <- NA_real_; prev_dt <- NA_real_

      seq_pos <- if (direction == 1) {
        (start + 1):ni
      } else {
        (start - 1):1
      }
      if (length(seq_pos) == 0) return(list(prob = probs, is_partial = partial))

      for (pos in seq_pos) {
        res <- score_triplet(prev_anchor, anchor, pos, prev_v, prev_w, prev_dt)
        probs[pos]   <- res$prob
        partial[pos] <- res$is_partial

        thr_pos <- if (isTRUE(res$is_partial)) threshold_partial else threshold_full
        if (!is.na(res$prob) && res$prob >= thr_pos) {
          prev_anchor <- anchor
          prev_v <- res$v
          prev_w <- res$w
          prev_dt <- res$dt
          anchor <- pos
        }

        if (abs(pos - anchor) > max_skip) {
          prev_anchor <- anchor
          prev_v <- NA_real_; prev_w <- NA_real_; prev_dt <- NA_real_
          anchor <- pos
        }
      }
      list(prob = probs, is_partial = partial)
    }

    ## Helper: per-fix threshold lookup given an is_partial mask.
    ## Vectorised; NA in mask -> NA in threshold.
    thr_for <- function(is_partial_vec) {
      ifelse(is_partial_vec, threshold_partial, threshold_full)
    }

    ## ---- apply scanning strategy ----
    ## Each strategy now also produces a `seq_ratio` (prob /
    ## regime-matched threshold) parallel to `seq_prob`.  Final flagging
    ## uses `seq_ratio < 1` rather than `seq_prob < threshold` so the
    ## two-threshold regime is correctly applied per fix.
    if (scan == "forward-backward") {
      fwd <- run_scan(1, +1)
      bwd <- run_scan(ni, -1)
      seq_prob <- pmax(fwd$prob, bwd$prob, na.rm = TRUE)
      ratio_fwd <- fwd$prob / thr_for(fwd$is_partial)
      ratio_bwd <- bwd$prob / thr_for(bwd$is_partial)
      seq_ratio <- pmax(ratio_fwd, ratio_bwd, na.rm = TRUE)

    } else if (scan == "greedy") {
      ## start from the highest-probability location
      sp <- sim_probs[idx]
      best <- which.max(sp)
      if (length(best) == 0) best <- 1
      right <- run_scan(best, +1)
      left  <- run_scan(best, -1)
      seq_prob <- pmax(right$prob, left$prob, na.rm = TRUE)
      ratio_right <- right$prob / thr_for(right$is_partial)
      ratio_left  <- left$prob  / thr_for(left$is_partial)
      seq_ratio <- pmax(ratio_right, ratio_left, na.rm = TRUE)
      ## the anchor itself gets its simultaneous prob (full-formula
      ## scope: mt_flag_outliers computes joint_prob with full autodiff).
      seq_prob[best]  <- sp[best]
      seq_ratio[best] <- sp[best] / threshold_full

    } else if (scan == "random") {
      ## launch scans from n_random anchors, take median assessment per point
      margin <- max(1, ni %/% 20)
      pool <- seq(margin, ni - margin)
      anchors <- sort(sample(pool, min(n_random, length(pool))))
      n_scans <- length(anchors) * 2
      prob_mat    <- matrix(NA_real_, nrow = ni, ncol = n_scans)
      partial_mat <- matrix(NA, nrow = ni, ncol = n_scans)

      for (k in seq_along(anchors)) {
        fwd_k <- run_scan(anchors[k], +1)
        bwd_k <- run_scan(anchors[k], -1)
        prob_mat   [, 2 * k - 1] <- fwd_k$prob
        partial_mat[, 2 * k - 1] <- fwd_k$is_partial
        prob_mat   [, 2 * k]     <- bwd_k$prob
        partial_mat[, 2 * k]     <- bwd_k$is_partial
      }
      thr_mat   <- ifelse(partial_mat, threshold_partial, threshold_full)
      ratio_mat <- prob_mat / thr_mat
      if (requireNamespace("matrixStats", quietly = TRUE)) {
        seq_prob  <- matrixStats::rowMedians(prob_mat,  na.rm = TRUE)
        seq_ratio <- matrixStats::rowMedians(ratio_mat, na.rm = TRUE)
      } else {
        seq_prob  <- apply(prob_mat,  1, stats::median, na.rm = TRUE)
        seq_ratio <- apply(ratio_mat, 1, stats::median, na.rm = TRUE)
      }
      seq_prob [is.nan(seq_prob)]  <- NA_real_
      seq_ratio[is.nan(seq_ratio)] <- NA_real_
    }

    ## ---- flag based on regime-matched threshold ----
    ## A fix passes if its best-evidence ratio (prob / matching threshold)
    ## >= 1.  Equivalent to the old `seq_prob < threshold` rule under
    ## the single-threshold regime when threshold_partial == threshold_full,
    ## so byte-identical on single-state tracks.
    is_out <- !is.na(seq_ratio) & seq_ratio < 1
    is_out[1] <- FALSE
    is_out[ni] <- FALSE

    all_outlier[idx] <- is_out
    all_prob[idx] <- seq_prob

    ## ---- anchor-corruption diagnostic --------------------------------
    ## A single corrupt anchor makes every legitimate fix look anomalous
    ## relative to it; below max_skip the failure mode is silent.  When
    ## the per-track flag rate exceeds the documented threshold (default
    ## 30%), warn the user that the scan is more likely reflecting an
    ## anchor problem than data contamination.
    if (!is.null(anchor_corruption_threshold)) {
      track_flag_rate <- sum(is_out) / ni
      if (track_flag_rate > anchor_corruption_threshold) {
        rlang::warn(
          sprintf(paste0(
            "Sequential scan flagged %.0f%% of fixes on track '%s' (%d of %d) ",
            "-- above the %.0f%% anchor-corruption threshold.  This is the ",
            "signature of a corrupt anchor: every legitimate fix looks ",
            "anomalous relative to a bad starting point.  Recovery: try ",
            "`scan = \"forward-backward\"` (default; uses both endpoints), ",
            "supply a clean `reference =` object, or use `mt_clean_track()` ",
            "which is anchor-free.  See `?mt_sequential_outliers` section ",
            "\"Anchor-corruption diagnostic\"."),
            100 * track_flag_rate, uid, sum(is_out), ni,
            100 * anchor_corruption_threshold),
          class = "move2utils_mt_sequential_outliers_anchor_corruption")
      }
    }
  }

  x$is_outlier <- all_outlier
  x$seq_joint_prob <- all_prob
  ## restore the caller's original CRS (flag columns are CRS-invariant)
  if (sf::st_crs(x) != orig_crs) x <- sf::st_transform(x, orig_crs)

  n_out <- sum(all_outlier)
  message(sprintf(
    "Sequential scan: flagged %d outliers (%.2f%% of %d locations).",
    n_out, 100 * n_out / nrow(x), nrow(x)
  ))

  if (plot && n_out > 0) {
    coords <- sf::st_coordinates(x)
    graphics::plot(coords, type = "n", asp = 1,
                   xlab = "Longitude", ylab = "Latitude",
                   main = paste("Sequential outlier detection:", scan))
    graphics::points(coords[!all_outlier, , drop = FALSE],
                     col = "grey60", pch = 19, cex = 0.3)
    graphics::points(coords[all_outlier, , drop = FALSE],
                     col = "red", pch = 4, cex = 0.8, lwd = 1.5)
  }

  x
}

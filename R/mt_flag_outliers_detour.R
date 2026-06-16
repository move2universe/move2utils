#' Flag outliers using path-vs-displacement detour ratio
#'
#' Time-insensitive, scale-invariant point-level outlier detector.
#' For each interior fix \eqn{N} and window radius \eqn{k}, computes
#' \deqn{ratio_k(N) = \frac{path\_length(N-k, \dots, N+k)}{displacement(N-k, N+k)}.}
#' Legitimate animal movement has \eqn{ratio_k \approx 1} when motion is
#' roughly straight and bounded above (typically \eqn{< 2}) even with
#' hard turns. An out-and-back GPS spike has displacement near zero but
#' path length preserved, so the ratio explodes -- regardless of the
#' time elapsed across the window.
#'
#' @details
#' **Why a separate primitive.** This metric is the topological
#' complement to the Brownian-bridge perpendicular residual produced
#' by \code{\link{mt_flag_outliers_bridge}}. Bridge residual is in
#' metres and scales with bridge length and Brownian variance; detour
#' is dimensionless and uses only successive coordinates. At sparse
#' sampling rates (e.g. 1-h GPS) the bridge \eqn{\sigma} broadens and
#' loses sensitivity to single-fix spikes whose implied step speed is
#' below physiological caps; detour stays robust because it only
#' compares spatial path to displacement. The two primitives flag
#' overlapping but distinct sets of fixes and combine well via the
#' conjunction rule in \code{\link{mt_clean_track}}.
#'
#' **Time-insensitivity.** The detector uses only the spatial sequence
#' of coordinates. It will produce the same flag whether your sampling
#' is uniform 1-h, irregular bursts, or post-resampled. This is
#' intentional and a strength when the kinematic primitives lose
#' precision under irregular sampling.
#'
#' **Window radius \code{k}.** \code{k = 1} catches single-fix spikes
#' (the dominant pattern at the colony in real GPS data). Larger
#' \code{k} catches block-shaped errors (2+ consecutive bad fixes
#' through which the bridge primitive can pass cleanly). Pass an
#' integer (single \code{k}) or a vector to compute the per-fix
#' \emph{maximum} ratio across a range of window sizes; in the latter
#' case the per-\code{k} threshold scaling matters because legitimate
#' movement at large \code{k} can produce ratio > 1.5 just from
#' mild detours, so prefer single \code{k = 1} for flagging and
#' larger \code{k} as diagnostic.
#'
#' **Leg gate \code{min_leg}.** Setting \code{min_leg > 0} requires
#' both incident step lengths to exceed \code{min_leg} metres before
#' a fix can be flagged. Prevents the detector from firing on
#' small-displacement noise wiggles where the ratio can be large but
#' the absolute path is tiny. When the detector is used inside
#' \code{\link{mt_clean_track}}, the conjunction with another
#' primitive plays the same gating role, so \code{min_leg = 0} is
#' usually fine in that context.
#'
#' @param x A \code{move2} object. Single- or multi-track. CRS:
#'   lon/lat is auto-handled via Haversine; otherwise coordinates
#'   are assumed to be in a projected CRS (metres).
#' @param k Integer (single) or integer vector. Window radius in
#'   fixes. Default \code{1L}. When a vector is supplied, the
#'   per-fix score is the maximum ratio across the supplied
#'   \code{k} values.
#' @param threshold Numeric. Detour ratio above which a fix is
#'   flagged.  Used only when \code{threshold_type = "fixed"}
#'   (the default).  Default \code{5}.  HEURISTIC --
#'   comfortably above the legitimate-flight ceiling for any
#'   sampling rate; plausible range 3--10.  Ignored when
#'   \code{threshold_type = "auto"}.
#' @param threshold_type Character, one of \code{"fixed"} (default,
#'   legacy) or \code{"auto"}.  In \code{"fixed"} mode, the
#'   user-supplied scalar \code{threshold} is used.  In
#'   \code{"auto"} mode the function computes the entropy valley
#'   in the empirical \eqn{-\log(\rho)} distribution -- outliers
#'   sit in the lower tail since they have \eqn{\rho \gg 1} and
#'   therefore \eqn{-\log(\rho) \ll 0} -- and uses the back-
#'   converted ratio (\eqn{\exp(-\text{break})}) as the threshold.
#'   Same threshold-detection machinery as the bridge / speed-cap /
#'   probability primitives' entropy paths
#'   (\code{\link{.entropy_threshold_lower}}, sweep-validated
#'   density-ratio of 0.3 unified package-wide 2026-05-09).  When
#'   no entropy valley is found, no fixes are flagged
#'   (entropy-detector's safe-on-clean contract).
#' @param min_leg Numeric, non-negative metres. Default \code{0}.
#'   When \code{> 0}, only flag fixes where both incident step
#'   lengths exceed \code{min_leg}. Useful for standalone use;
#'   leave at \code{0} when combined with another primitive via
#'   the conjunction rule.
#' @param pool_by Optional character vector of length 1 or 2 naming
#'   column(s) in \code{mt_track_data(x)}.  Length 1: single column
#'   used as both fit set and operating unit (e.g.
#'   \code{"individual_id"}).  Length 2: \code{c(outer, inner)}
#'   where \code{outer} names the fit-source column (the union of
#'   its events supplies the entropy-break distribution) and
#'   \code{inner} names the operating unit (within which pool-added
#'   flags are unioned).  Length 2 requires strict nesting: every
#'   distinct \code{inner} value must map to exactly one
#'   \code{outer} value.  Length \eqn{> 2} is rejected -- pool_by
#'   has exactly two semantic roles, and deeper hierarchies would
#'   require hierarchical threshold estimation, which this primitive
#'   does not perform.  Pool-fit flags are unioned into per-track
#'   flags -- the pool path can only add flags, never remove ones
#'   the per-track pass caught.  \code{NULL} (default) preserves
#'   per-track behaviour byte-identically.  NA values in the named
#'   column(s) cause those tracks to fall back to per-track
#'   processing with a warning.  Only relevant for
#'   \code{threshold_type = "auto"}; \code{"fixed"} ignores
#'   \code{pool_by} (the user-supplied scalar is already the same
#'   across tracks).  See \code{?mt_clean_track} for the full
#'   semantics of the orchestrator's post-cascade pool sweep.
#' @param plot Logical. Diagnostic map of flagged fixes. Default
#'   \code{TRUE}.
#' @param remove Logical. If \code{TRUE}, return only kept rows.
#'   Default \code{FALSE} (return all rows with flag columns
#'   attached).
#' @param silent Logical. Suppress narration. Default \code{FALSE}.
#'
#' @return A \code{move2} object with added columns:
#'   \describe{
#'     \item{\code{is_outlier}}{Logical. TRUE where flagged.}
#'     \item{\code{flagged_by_detour}}{Logical. Same as
#'       \code{is_outlier} for this primitive; named for parity with
#'       the other primitives' output schema.}
#'     \item{\code{detour_ratio}}{Numeric. Per-fix maximum ratio
#'       across the supplied \code{k} values; \code{NA} at boundary
#'       fixes where the window cannot fit.}
#'     \item{\code{loglr_detour}}{Signed log-likelihood-ratio of outlier
#'       vs not, in nat units (surprisal \code{log(detour_ratio)} minus
#'       the detector's own flag boundary); \code{> 0} where flagged,
#'       \code{NA} where there is no detour geometry (ratio \eqn{\le} 1).
#'       See \code{DESIGN_evidence_accumulation.md}.}
#'   }
#'
#' @seealso \code{\link{mt_flag_outliers_bridge}},
#'   \code{\link{mt_flag_outliers}},
#'   \code{\link{mt_flag_speed_cap}},
#'   \code{\link{mt_clean_track}}.
#'
#' @examples
#' \dontrun{
#' library(move2)
#' x <- movebank_download_study(study_id = 123)
#' x <- mt_filter_gps_quality(x)
#' out <- mt_flag_outliers_detour(x, k = 1, threshold = 5,
#'                                 min_leg = 5000)
#' table(out$is_outlier)
#' }
#'
#' @importFrom move2 mt_track_id
#' @importFrom sf st_coordinates st_is_longlat
#' @export
mt_flag_outliers_detour <- function(x,
                                     k = 1L,
                                     threshold = 5,
                                     threshold_type = c("fixed", "auto"),
                                     min_leg = 0,
                                     pool_by = NULL,
                                     plot = TRUE,
                                     remove = FALSE,
                                     silent = FALSE) {
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  if (!is.numeric(k) || any(k < 1) || any(k != as.integer(k))) {
    rlang::abort("`k` must be a positive integer or vector of positive integers.",
                 class = "move2utils_mt_flag_outliers_detour_bad_k")
  }
  k <- as.integer(k)
  threshold_type <- match.arg(threshold_type)
  ## `threshold` is only validated as a numeric > 1 in the "fixed"
  ## path; in the "auto" path the value is replaced by an entropy-
  ## valley break detected on -log(ratio), and the input scalar is
  ## ignored.  Same convention used by mt_flag_speed_cap's
  ## threshold_type = "auto" / "hard" branches.
  if (threshold_type == "fixed") {
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        is.na(threshold) || threshold <= 1) {
      rlang::abort("`threshold` must be a scalar > 1.",
                   class = "move2utils_mt_flag_outliers_detour_bad_threshold")
    }
  }
  if (!is.numeric(min_leg) || length(min_leg) != 1L ||
      is.na(min_leg) || min_leg < 0) {
    rlang::abort("`min_leg` must be a non-negative scalar (metres).",
                 class = "move2utils_mt_flag_outliers_detour_bad_min_leg")
  }
  if (!is.null(pool_by)) {
    ## `.resolve_pool_groups` validates shape, columns, and (for
    ## length-2 input) the strict-nesting requirement.  Call early so
    ## user-facing errors fire on bad inputs regardless of option
    ## ordering (e.g. threshold_type = "fixed" later short-circuits
    ## the pool path).
    invisible(.resolve_pool_groups(x, pool_by, silent = TRUE))
  }
  was_longlat <- isTRUE(sf::st_is_longlat(x))

  ## ---- common plumbing: multi-track dispatch + extraction + hygiene
  ## + .detour_fn_core + lift.  Skip projection: detour intentionally
  ## stays in input CRS (Haversine on lon/lat, Euclidean on projected)
  ## to avoid AEQD projection cost on multi-million-fix tracks. ----
  out <- .clean_track_dispatch(
    x,
    fn_core             = .detour_fn_core,
    fn_core_args        = list(was_longlat   = was_longlat,
                                k             = k,
                                threshold     = threshold,
                                threshold_type = threshold_type,
                                min_leg       = min_leg,
                                silent        = silent),
    per_track_args      = list(),
    need_time           = FALSE,
    hygiene_strict_time = FALSE,
    project_longlat     = FALSE,
    n_min               = 2L * max(k) + 1L,
    n_min_severity      = "say",
    primitive_label     = "detour",
    silent              = silent)

  ## ---- pool_by union: refit threshold per group and union flags ----
  ## Only runs when pool_by is set AND threshold_type = "auto".  The
  ## "fixed" path uses a user-supplied scalar that is already the
  ## same across tracks; pooling does nothing there.
  if (!is.null(pool_by) && threshold_type == "auto") {
    out <- .detour_pool_union(out, x, pool_by, k = k,
                                min_leg = min_leg,
                                was_longlat = was_longlat,
                                silent = silent)
  }

  ## log-LR emission (evidence currency; additive, is_outlier untouched)
  out <- .attach_detour_loglr(out)

  if (plot)   .plot_detour(out)
  if (remove) out <- out[!out$is_outlier, ]
  out
}

## Detour log-LR.  detour_ratio = path / displacement (>= 1, higher =
## more out-and-back = more outlier; the detector thresholds
## -log(ratio)).  The surprisal is log(ratio) nats -- 0 at no detour
## (ratio 1), rising as the path inflates over the chord.  Ratios <= 1 or
## non-finite have no detour geometry -> abstain (NA).
## @keywords internal
.attach_detour_loglr <- function(obj) {
  r <- obj$detour_ratio
  if (is.null(r)) return(obj)
  neglogp <- log(r)
  neglogp[!is.finite(neglogp) | !is.finite(r) | r <= 1] <- NA_real_
  obj$loglr_detour <- .loglr_grouped(neglogp, obj$is_outlier,
                                     move2::mt_track_id(obj))
  obj
}


## Pool-union closure factory for mt_flag_outliers_detour.
##
## For each pool group, fits one entropy threshold from the union
## of all group tracks' detour_ratio values, then re-flags each
## track using the pool threshold (with the same legs_ok gate the
## per-track pass applied).  Unions into is_outlier +
## flagged_by_detour.
##
## @keywords internal
.detour_pool_union <- function(out, x, pool_by, k, min_leg,
                                was_longlat, silent = FALSE) {

  pool_step <- function(fit_idx, apply_idx) {
    ## Fit pool threshold from the OUTER group's union of detour ratios.
    ratio_fit <- out$detour_ratio[fit_idx]
    valid <- !is.na(ratio_fit) & is.finite(ratio_fit) & ratio_fit > 1
    if (sum(valid) < 10L) return(rep(FALSE, length(apply_idx)))
    nl <- rep(NA_real_, length(ratio_fit))
    nl[valid] <- -log(ratio_fit[valid])
    ## Defer to .entropy_threshold_lower's leaf formal -- single source
    ## of truth for the package-wide entropy default.  Pre-2026-05-25
    ## the literal 0.3 was forwarded explicitly here, shadowing the
    ## leaf.  See audits/2026-05-25-parameter-propagation/findings.md §1.
    res <- .entropy_threshold_lower(nl)
    if (is.na(res$break_value)) return(rep(FALSE, length(apply_idx)))
    pool_threshold <- exp(-res$break_value)

    ## Apply to INNER group: each apply_idx event's detour ratio.
    ratio_g <- out$detour_ratio[apply_idx]

    ## legs_ok gate -- if min_leg > 0, recompute per-track step
    ## lengths and require both incident steps >= min_leg, matching
    ## the gate inside .detour_fn_core.  When min_leg = 0, no gate.
    legs_ok <- if (min_leg > 0) {
      ids_g  <- as.character(move2::mt_track_id(out))[apply_idx]
      out_ok <- rep(TRUE, length(apply_idx))
      for (tid in unique(ids_g)) {
        rel <- which(ids_g == tid)
        cc_t <- sf::st_coordinates(out[apply_idx[rel], ])
        sd_full <- .step_lengths_from_cc(cc_t, was_longlat)
        step_in  <- c(NA_real_, sd_full[-length(sd_full)])
        step_out <- sd_full
        out_ok[rel] <- !is.na(step_in) & !is.na(step_out) &
                        step_in >= min_leg & step_out >= min_leg
      }
      out_ok
    } else {
      rep(TRUE, length(apply_idx))
    }

    flagged <- !is.na(ratio_g) & ratio_g > pool_threshold & legs_ok
    if (!silent) {
      n_new <- sum(flagged & !out$is_outlier[apply_idx], na.rm = TRUE)
      message(sprintf(
        "  pool_by[%s]: threshold = %.3g, %d new fix(es) flagged across %d event(s).",
        paste(pool_by, collapse = ","), pool_threshold, n_new, length(apply_idx)))
    }
    flagged
  }

  .apply_pool_union(out, x, pool_by, pool_step,
                    flag_cols = c("is_outlier", "flagged_by_detour"),
                    silent = silent)
}


## Raw-matrix entry point used by mt_clean_track to skip per-iter
## sf-class slicing.  Operates on (cc, active_idx, was_longlat, k, ...)
## and returns active-indexed result vectors.
##
## Detour intentionally avoids AEQD projection: path and displacement
## are harmonised on the same metric (great-circle for lon/lat input,
## Euclidean for projected input).  This avoids the O(n) projection
## cost on multi-million-fix tracks where the avoided projection is
## the single dominant work item.
##
## @keywords internal
.detour_fn_core <- function(cc, active_idx, was_longlat,
                             k              = 1L,
                             threshold      = 5,
                             threshold_type = "fixed",
                             min_leg        = 0,
                             silent         = TRUE) {

  ## Contract: active_idx sorted ascending with no duplicates.
  stopifnot(!is.unsorted(active_idx), !anyDuplicated(active_idx))

  say   <- function(...) if (!silent) message(...)
  k     <- as.integer(k)
  k_max <- max(k)
  n_a   <- length(active_idx)

  if (n_a < 2L * k_max + 1L) {
    say("Track too short for detour computation; returning unflagged.")
    return(list(
      is_outlier            = rep(FALSE, n_a),
      flagged_by_detour     = rep(FALSE, n_a),
      detour_ratio          = rep(NA_real_, n_a),
      detour_threshold_used = NA_real_
    ))
  }

  cc_a <- cc[active_idx, , drop = FALSE]

  ## Step distances on the active subset.  The last element is NA per
  ## the package convention; replace with 0 for the cumulative-path
  ## running sum.
  step_dist <- .step_lengths_from_cc(cc_a, was_longlat)
  step_dist[is.na(step_dist)] <- 0
  cum_path <- c(0, cumsum(step_dist[-n_a]))

  ## Per-k ratio computation, vectorised.
  ratios <- matrix(NA_real_, nrow = n_a, ncol = length(k))
  for (j in seq_along(k)) {
    kj  <- k[j]
    idx <- (kj + 1L):(n_a - kj)
    path_k <- cum_path[idx + kj] - cum_path[idx - kj]
    if (was_longlat) {
      disp_k <- .haversine_pair(cc_a[idx - kj, 1], cc_a[idx - kj, 2],
                                  cc_a[idx + kj, 1], cc_a[idx + kj, 2])
    } else {
      dx <- cc_a[idx + kj, 1] - cc_a[idx - kj, 1]
      dy <- cc_a[idx + kj, 2] - cc_a[idx - kj, 2]
      disp_k <- sqrt(dx ^ 2 + dy ^ 2)
    }
    ratios[idx, j] <- path_k / pmax(disp_k, 1)  # avoid /0
  }
  ratio <- if (length(k) == 1L) ratios[, 1L] else
    suppressWarnings(apply(ratios, 1, max, na.rm = TRUE))
  ratio[is.infinite(ratio)] <- NA_real_

  ## Leg gate (only when min_leg > 0).
  if (min_leg > 0) {
    sd_full  <- .step_lengths_from_cc(cc_a, was_longlat)
    step_in  <- c(NA_real_, sd_full[-n_a])
    step_out <- sd_full
    legs_ok  <- !is.na(step_in) & !is.na(step_out) &
                step_in >= min_leg & step_out >= min_leg
  } else {
    legs_ok <- rep(TRUE, n_a)
  }

  ## ---- resolve threshold ----
  ## "fixed" (default): user-supplied scalar.
  ## "auto": entropy valley in -log(ratio); convert back to ratio
  ## space.  No valley -> Inf threshold (no flags; entropy detector's
  ## safe-on-clean contract).  Uses the package-wide entropy default
  ## (sweep-validated 2026-05-06; defaults to
  ## .entropy_threshold_lower's leaf formal).
  if (threshold_type == "auto") {
    valid_ratio <- !is.na(ratio) & is.finite(ratio) & ratio > 1
    if (sum(valid_ratio) >= 10L) {
      neg_log_ratio <- rep(NA_real_, length(ratio))
      neg_log_ratio[valid_ratio] <- -log(ratio[valid_ratio])
      ## Leaf-formal default (single source of truth); see
      ## audits/2026-05-25-parameter-propagation/findings.md §1.
      res_thr <- .entropy_threshold_lower(neg_log_ratio)
      threshold_resolved <- if (!is.na(res_thr$break_value)) {
        exp(-res_thr$break_value)
      } else {
        Inf
      }
    } else {
      threshold_resolved <- Inf
    }
  } else {
    threshold_resolved <- threshold
  }

  flagged <- !is.na(ratio) & ratio > threshold_resolved & legs_ok

  thr_msg <- if (threshold_type == "auto") {
    if (is.finite(threshold_resolved)) {
      sprintf("auto -> %.3g", threshold_resolved)
    } else {
      "auto -> Inf (no entropy valley found; safe-on-clean default)"
    }
  } else {
    sprintf("%g", threshold_resolved)
  }
  say(sprintf("Detour primitive: flagged %d of %d fixes (k=%s, threshold=%s).",
              sum(flagged), n_a,
              if (length(k) == 1L) as.character(k) else
                paste0("{", paste(k, collapse = ","), "}"),
              thr_msg))

  list(
    is_outlier            = flagged,
    flagged_by_detour     = flagged,
    detour_ratio          = ratio,
    detour_threshold_used = threshold_resolved
  )
}


## Diagnostic plot for mt_flag_outliers_detour() output.
##
## @keywords internal
.plot_detour <- function(x) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(3, 3, 3, 1))

  cc <- sf::st_coordinates(x)
  graphics::plot(cc, type = "l", col = "grey70", lwd = 0.3,
                 asp = 1, xlab = "", ylab = "",
                 main = sprintf(
                   "mt_flag_outliers_detour -- %d flagged (%.2f%%)",
                   sum(x$is_outlier, na.rm = TRUE),
                   100 * mean(x$is_outlier, na.rm = TRUE)))
  flagged <- which(x$is_outlier)
  if (length(flagged)) {
    graphics::points(cc[flagged, , drop = FALSE],
                     col = grDevices::adjustcolor("red", 0.6),
                     pch = 20, cex = 0.6)
  }
  invisible(NULL)
}

#' Multi-scale persistence score for outlier flags
#'
#' Annotates flagged fixes with a \emph{persistence score} -- the
#' number of temporal scales at which each flag's local geometry is
#' independently anomalous.  The score is detector-agnostic: it
#' operates on any \code{move2} object that carries an \code{is_outlier}
#' column, regardless of which primitive flagged the fixes
#' (\code{\link{mt_clean_track}}, \code{\link{mt_flag_outliers_bridge}},
#' \code{\link{mt_flag_outliers_detour}}, \code{\link{mt_flag_outliers}},
#' \code{\link{mt_flag_speed_cap}}, etc.).
#'
#' @details
#' \strong{Mechanism.}  For each flagged fix \eqn{i} and each scale
#' \eqn{k}, the function computes a "scale-\eqn{k} view" of \eqn{i}'s
#' local geometry: the step lengths from \eqn{x_{i-k}} to \eqn{x_i}
#' and from \eqn{x_i} to \eqn{x_{i+k}}, and the turn angle at \eqn{i}
#' between those two long arms.  The reference distribution at scale
#' \eqn{k} is built from the same scale-\eqn{k} geometry computed at
#' every interior fix of the track.  A 2-D histogram of
#' \eqn{(\log\text{step}_k, \text{turn}_k)} is the joint density; the
#' flagged fix's bin density gives its scale-\eqn{k} probability; a
#' gap-on-(-log probability) threshold determines whether the fix is
#' anomalous \emph{at scale \eqn{k}}.
#'
#' Crucially, the scale-\eqn{k} view is computed for \emph{every}
#' fix, not just for fixes that would be retained in a thinned-track
#' grid.  This dissolves the index-parity issue of conventional
#' multi-scale voting (where some fixes can never be evaluated at
#' coarser scales because the thinning grid skips them).
#'
#' \strong{Persistence score.}  A flagged fix's persistence score is
#' \deqn{p(i) = 1 + \sum_{k \in \text{scales}} \mathbb{1}[i \text{ flagged at scale } k],}
#' where the constant 1 accounts for scale 1 (the original detector's
#' flag).  With the default \code{scales = c(2, 4, 8)} the score
#' takes integer values from 1 (flagged only at native resolution)
#' to 4 (flagged at every scale tested).  Higher scores indicate
#' anomalies whose geometric extent survives temporal coarsening.
#'
#' \strong{Class-aware filtering recommendation.}  Empirical work on
#' the synthetic CPF ground truth (see vignette) shows persistence
#' is class-conditionally informative when the input is cascade
#' output (\code{\link{mt_clean_track}}):
#'
#' \itemize{
#'   \item \code{geometric_spike} class: empirically class-pure on
#'     the synthetic; persistence has nothing to filter.
#'   \item \code{state_anomaly} and \code{consensus} classes: TPs
#'     persist 80--95\% at \code{p >= 3}, FPs only 44--57\%.  A
#'     filter at \code{persistence_count >= 3} substantively
#'     improves precision on these classes.
#'   \item \code{kinematic_confluence} class: small sample on the
#'     synthetic; signal direction is uncertain.
#' }
#'
#' Recommendation: use the score as a \emph{confidence column}, not
#' an automatic filter.  Where a filter is desired, gate it on
#' \code{error_class} (when the input came from
#' \code{\link{mt_clean_track}}).
#'
#' \strong{When to expect persistence to be informative.}  The
#' persistence score is most useful when the underlying error type
#' has \emph{geometric extent that survives coarsening}: isolated
#' single-fix spikes (whose 4-step or 8-step neighbourhood is still
#' dominated by the spike), and block-boundary fixes (where the
#' coarse-scale step crosses the spoof boundary).  It is least
#' useful for halo-style outliers (repeated wandering returns) whose
#' anomaly averages out over wider windows.
#'
#' @param x A \code{move2} object with an \code{is_outlier} logical
#'   column (single- or multi-track).
#' @param scales Integer vector of validation scales \eqn{k}.  Each
#'   \eqn{k} must satisfy \eqn{2 \le k \le \lfloor (n-1)/2 \rfloor}
#'   on a track with \eqn{n} fixes (boundary fixes near the track
#'   ends cannot be evaluated at scale \eqn{k}).  Default
#'   \code{c(2, 4, 8)}.
#' @param threshold Numeric.  Gap-threshold parameter passed to the
#'   internal \code{.gap_threshold_lower} for per-scale flag
#'   detection on \eqn{-\log(\text{prob})}.  Default \code{NULL}, which
#'   defers to the leaf default of \code{.gap_threshold_lower} (3).
#' @param n_breaks Integer.  Number of bins per axis for the 2-D
#'   per-scale joint-probability histogram.  Default \code{20L}.
#' @param silent Logical.  If \code{TRUE} suppress per-track
#'   summary messages.  Default \code{FALSE}.
#'
#' @return The input \code{x} with added columns:
#'   \describe{
#'     \item{\code{persistence_count}}{Integer.  The persistence
#'       score (1 = flagged only at native resolution; up to
#'       \code{1 + length(scales)} = flagged at every validation
#'       scale).  \code{NA_integer_} for fixes that were not flagged
#'       in \code{is_outlier} (the score is defined only for
#'       candidates).}
#'     \item{\code{persistence_at_scale_<k>}}{One logical column per
#'       scale \eqn{k} in \code{scales}, indicating whether the
#'       flagged fix's scale-\eqn{k} geometry was anomalous against
#'       the scale-\eqn{k} reference distribution.  \code{NA} for
#'       non-flagged fixes.}
#'   }
#'
#' \code{is_outlier} itself is \emph{not} modified.  The function is
#' a pure annotator -- users decide whether and how to filter.
#'
#' @examples
#' \dontrun{
#' ## Cascade output annotated with persistence scores
#' clean <- mt_clean_track(track, mass = 5, mode = "flying", remove = FALSE)
#' annotated <- mt_persistence_score(clean)
#'
#' ## Class-aware filter (recommended pattern):
#' is_low_confidence <-
#'   annotated$is_outlier &
#'   annotated$error_class %in% c("state_anomaly", "consensus") &
#'   annotated$persistence_count < 3
#' annotated$is_outlier[is_low_confidence] <- FALSE
#' }
#'
#' @seealso \code{\link{mt_clean_track}} for the cascade orchestrator
#'   whose \code{error_class} column gates the filtering rule;
#'   \code{vignette("OUTLIER_5_persistence_score", package = "move2utils")} for
#'   the empirical class-conditional analysis on synthetic CPF data.
#'
#' @importFrom move2 mt_track_id mt_aeqd_crs
#' @importFrom sf st_coordinates st_is_longlat st_transform
#' @export
mt_persistence_score <- function(x,
                                   scales    = c(2L, 4L, 8L),
                                   threshold = NULL,
                                   n_breaks  = 20L,
                                   silent    = FALSE) {
  ## ---- input validation ------------------------------------------
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  if (!"is_outlier" %in% names(x)) {
    rlang::abort(paste0(
      "`x` must carry an `is_outlier` column from a flagger ",
      "(e.g. mt_clean_track, mt_flag_outliers_bridge)."),
      class = "move2utils_mt_persistence_score_missing_is_outlier")
  }
  if (!is.logical(x$is_outlier)) {
    rlang::abort("`is_outlier` must be a logical column.",
                 class = "move2utils_mt_persistence_score_is_outlier_not_logical")
  }
  scales <- as.integer(scales)
  if (length(scales) < 1L || any(is.na(scales)) || any(scales < 2L)) {
    rlang::abort(paste0(
      "`scales` must be a non-empty integer vector with all ",
      "values >= 2."),
      class = "move2utils_mt_persistence_score_bad_scales")
  }
  ## NULL defers to .persistence_fn_core's leaf formal (single source
  ## of truth for the per-scale gap-threshold default).  Non-NULL
  ## scalar overrides.  See
  ## audits/2026-05-25-parameter-propagation/findings.md §1.10.
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        is.na(threshold) || threshold <= 0) {
      rlang::abort("`threshold` must be a positive scalar or NULL.",
                   class = "move2utils_mt_persistence_score_bad_threshold")
    }
  }
  n_breaks <- as.integer(n_breaks)
  if (is.na(n_breaks) || n_breaks < 4L) {
    rlang::abort("`n_breaks` must be an integer >= 4.",
                 class = "move2utils_mt_persistence_score_bad_n_breaks")
  }

  ## ---- common plumbing: multi-track dispatch + AEQD projection
  ## (deterministic histogram bins) + extraction + .persistence_fn_core
  ## + lift.  The candidate set is x$is_outlier, passed as a per-track
  ## sliced logical vector.  NULL threshold is omitted from
  ## fn_core_args so the leaf formal is the single source of truth. ----
  fn_core_args <- list(scales = scales, n_breaks = n_breaks, silent = silent)
  if (!is.null(threshold)) fn_core_args$threshold <- threshold
  out <- .clean_track_dispatch(
    x,
    fn_core             = .persistence_fn_core,
    fn_core_args        = fn_core_args,
    per_track_args      = list(candidate_mask = x$is_outlier),
    need_time           = FALSE,
    hygiene_strict_time = FALSE,
    project_longlat     = TRUE,
    n_min               = 1L,
    n_min_severity      = "say",
    primitive_label     = "persistence",
    silent              = silent)

  out
}


## Raw-matrix entry point for multi-scale persistence annotation.
##
## Persistence is a geometric, time-insensitive consistency check on
## a candidate set (`candidate_mask`).  For each scale k in `scales`,
## the per-fix (step_in_k, step_out_k, turn_k) geometry is computed
## and a gap-threshold detector flags the candidates whose coarser-
## scale geometry is inconsistent with the body of the track.  Score
## = 1 (scale-1, by definition) + count of scales at which the
## candidate persists.
##
## Stays on projected coords (dispatcher's project_longlat=TRUE
## upstream guarantees a metric CRS) so the histogram bin breaks are
## CRS-deterministic.
##
## @keywords internal
.persistence_fn_core <- function(cc, active_idx, candidate_mask,
                                  scales    = c(2L, 4L, 8L),
                                  threshold = NULL,
                                  n_breaks  = 20L,
                                  silent    = TRUE) {
  ## NULL threshold defers to .gap_threshold_lower's leaf formal (via
  ## .flag_at_scale).  Pattern A; see
  ## audits/2026-05-25-parameter-propagation/findings.md §1.10.

  ## Contract: active_idx sorted ascending with no duplicates.
  stopifnot(!is.unsorted(active_idx), !anyDuplicated(active_idx))

  say   <- function(...) if (!silent) message(...)
  n_a   <- length(active_idx)
  out_cols <- paste0("persistence_at_scale_", scales)

  ## Allocate active-indexed return vectors.  Convention: NA outside
  ## the candidate set; integers / TRUE-FALSE at candidate positions.
  persistence_count <- rep(NA_integer_, n_a)
  flag_mat_active   <- vector("list", length(scales))
  names(flag_mat_active) <- out_cols
  for (cn in out_cols) flag_mat_active[[cn]] <- rep(NA, n_a)

  if (n_a == 0L) {
    return(c(list(persistence_count = persistence_count), flag_mat_active))
  }

  ## Restrict to candidates within the active subset.
  cand_local <- which(candidate_mask[active_idx])
  if (length(cand_local) == 0L) {
    say("No flagged fixes -- persistence score is defined only on candidates.")
    return(c(list(persistence_count = persistence_count), flag_mat_active))
  }

  ## Per-scale geometry computed on the active-subset coords; flags
  ## evaluated at the local candidate positions.
  cc_a <- cc[active_idx, , drop = FALSE]
  flag_mat <- matrix(FALSE, nrow = n_a, ncol = length(scales))
  colnames(flag_mat) <- as.character(scales)
  for (j in seq_along(scales)) {
    k <- scales[j]
    geom <- .scale_k_geometry(cc_a, k)
    flag_mat[, j] <- .flag_at_scale(geom, cand_local,
                                      threshold = threshold,
                                      n_breaks  = n_breaks)
  }

  persistence_count[cand_local] <-
    1L + as.integer(rowSums(flag_mat[cand_local, , drop = FALSE]))
  for (j in seq_along(scales)) {
    cn <- out_cols[j]
    flag_mat_active[[cn]][cand_local] <- flag_mat[cand_local, j]
  }

  pcounts <- persistence_count[cand_local]
  say(sprintf(
    "Persistence: %d flag(s) annotated; mean score = %.2f, max = %d.",
    length(cand_local), mean(pcounts), max(pcounts)))

  c(list(persistence_count = persistence_count), flag_mat_active)
}


## ---------------------------------------------------------------------
## Internal helpers
## ---------------------------------------------------------------------

## Compute (step_in_k, step_out_k, turn_k) at every interior fix of
## a single-track coordinate matrix where i - k >= 1 and i + k <= n.
## Returns a list of three numeric vectors of length n; NA where i is
## too close to a boundary.
##
## Coordinates are assumed to be in a planar projection (metres).
## mt_persistence_score auto-projects lon/lat input to AEQD before
## calling this helper, so the lon/lat dispatch that earlier versions
## carried is no longer needed.
##
## @keywords internal
.scale_k_geometry <- function(cc, k) {
  n <- nrow(cc)
  step_in  <- rep(NA_real_, n)
  step_out <- rep(NA_real_, n)
  turn     <- rep(NA_real_, n)
  if (n <= 2L * k) return(list(step_in = step_in,
                                step_out = step_out, turn = turn))

  i <- (k + 1L):(n - k)
  dx_in  <- cc[i,     1L] - cc[i - k, 1L]
  dy_in  <- cc[i,     2L] - cc[i - k, 2L]
  dx_out <- cc[i + k, 1L] - cc[i,     1L]
  dy_out <- cc[i + k, 2L] - cc[i,     2L]
  step_in[i]  <- sqrt(dx_in^2  + dy_in^2)
  step_out[i] <- sqrt(dx_out^2 + dy_out^2)
  turn[i]     <- atan2(dx_in * dy_out - dy_in * dx_out,
                        dx_in * dx_out + dy_in * dy_out)

  list(step_in = step_in, step_out = step_out, turn = turn)
}

## Per-scale flag rule.  Reference distribution: ALL fixes' scale-k
## (step_out, turn) values (including candidates -- self-contamination
## is mild because candidate count is small relative to track size,
## and excluding them complicates the API for marginal benefit).
## Threshold: gap-on-(-log prob).  Returns a logical vector of length
## nrow(geometry); TRUE only for indices in candidate_idx whose
## scale-k probability falls below the gap threshold.
##
## @keywords internal
.flag_at_scale <- function(geom, candidate_idx, threshold = NULL,
                              n_breaks = 20L) {
  step <- geom$step_out
  turn <- geom$turn
  ok <- is.finite(step) & is.finite(turn) & step > 0
  if (sum(ok) < 50L) {
    return(rep(FALSE, length(step)))
  }

  log_step <- log(step[ok])
  turn_ok  <- turn[ok]
  br_step <- seq(min(log_step), max(log_step), length.out = n_breaks + 1L)
  br_turn <- seq(-pi, pi,                      length.out = n_breaks + 1L)
  h <- table(
    cut(log_step, br_step, include.lowest = TRUE),
    cut(turn_ok,  br_turn,  include.lowest = TRUE))
  dens <- as.matrix(h) / sum(h)
  dens[dens == 0] <- .Machine$double.xmin

  bin_step <- rep(NA_integer_, length(step))
  bin_turn <- rep(NA_integer_, length(turn))
  bin_step[ok] <- as.integer(cut(log_step, br_step, include.lowest = TRUE))
  bin_turn[ok] <- as.integer(cut(turn_ok,  br_turn,  include.lowest = TRUE))

  prob <- rep(NA_real_, length(step))
  ok2 <- !is.na(bin_step) & !is.na(bin_turn)
  prob[ok2] <- dens[cbind(bin_step[ok2], bin_turn[ok2])]

  nl <- -log(prob)
  ok3 <- is.finite(nl)
  if (sum(ok3) < 10L) return(rep(FALSE, length(step)))
  ## NULL threshold defers to .gap_threshold_lower's leaf formal.
  thr <- if (is.null(threshold))
           .gap_threshold_lower(nl[ok3])$break_value
         else
           .gap_threshold_lower(nl[ok3], threshold = threshold)$break_value
  if (is.na(thr)) return(rep(FALSE, length(step)))

  flagged <- !is.na(nl) & nl > thr
  out <- rep(FALSE, length(step))
  out[candidate_idx] <- flagged[candidate_idx]
  out
}

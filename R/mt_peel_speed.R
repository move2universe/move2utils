#' Iterative speed peel at a fixed physiological cap
#'
#' Remove fixes whose implied step speed exceeds a user-supplied cap,
#' iteratively recomputing step speeds on the survivors so that peeling
#' a single fix exposes its former neighbours' step speeds and the
#' process can continue into coherent multi-fix error clusters.  Unlike
#' \code{\link{mt_flag_speed_cap}}, which applies the cap once to the
#' original distribution, \code{mt_peel_speed} keeps peeling until no
#' surviving fix has \code{max(step_in, step_out) > v_max}.
#'
#' @details
#' This is the primitive that catches coherent multi-fix errors that
#' per-fix detectors cannot resolve.  A K02-style spoof of 175 fixes
#' at a fake location has small step speeds internally (the fake
#' trajectory is self-consistent) and only the jump-out / jump-back
#' boundaries violate a physiological cap.  Peeling those boundaries
#' exposes the next interior fix as a new boundary (pre-spoof neighbour
#' now connects to an interior spoof fix over the full jump distance),
#' which also violates the cap and gets peeled.  Iteration walks inward
#' from both ends of the spoof until the clusters meet and the
#' remaining track is clean.
#'
#' \code{v_max} must be supplied by the user from species biology --
#' typically the maximum sustained flight speed or a clear
#' physiologically-impossible cap.  A data-driven \code{v_max} is
#' intentionally not offered here; see \code{\link{mt_suggest_speed_cap}}
#' for a diagnostic helper that inspects the distribution and reports
#' candidate cuts without choosing for you.
#'
#' Distances are computed in the input CRS when projected, or in metres
#' via per-track local AEQD projection when the input is longitude /
#' latitude.  The \code{step_speed} returned is always in m/s.
#'
#' @param x A \code{move2} object (single-track or multi-track).
#' @param v_max Numeric scalar, the physiological speed cap in m/s.
#'   Must be positive; no default.
#' @param aux_scores Optional numeric vector of auxiliary per-fix
#'   outlier scores aligned to \code{nrow(x)}.  When supplied, the peel
#'   switches to \emph{asymmetric} mode: instead of flagging both
#'   endpoints of every offending edge (\code{step_speed > v_max}), the
#'   function flags only the higher-scoring endpoint per edge.  Higher
#'   score = more likely outlier.  Typical sources: \code{bridge_eta}
#'   from \code{\link{mt_flag_outliers_bridge}}; \code{detour_ratio}
#'   from \code{\link{mt_flag_outliers_detour}}; \code{-log(prob)} from
#'   \code{\link{mt_flag_outliers}}; or any user-defined combination.
#'   On a 1-fix spike with reliable auxiliary scores, asymmetric mode
#'   removes the spike and preserves both clean neighbours where the
#'   default symmetric mode would remove all three.  Default
#'   \code{NULL} (symmetric mode).  See \emph{Asymmetric peel} below.
#' @param max_iter Integer.  Hard safety cap on peel iterations.
#'   Default 1000; in practice convergence is reached in tens of
#'   iterations for typical error clusters.
#' @param remove Logical.  If \code{TRUE}, return only surviving fixes
#'   (\code{is_outlier == FALSE}).  Default \code{FALSE}.
#' @param silent Logical.  If \code{FALSE} (default) the function
#'   prints a brief summary (cap used, fixes peeled, iterations).
#'   Set \code{TRUE} to suppress.  Errors and warnings are always
#'   shown.
#'
#' @section Asymmetric peel:
#' The default symmetric peel labels every fix whose
#' \code{max(step_in, step_out) > v_max} as an outlier in the current
#' iteration.  On a 1-fix spike both incident edges violate, so the
#' spike fix \emph{plus} its two clean neighbours all carry an
#' edge-violation flag and are removed in the same pass.  Only the
#' spike is structurally implausible; the two neighbours are clean
#' fixes whose only crime is being adjacent to the spike.
#'
#' When \code{aux_scores} is supplied, the function instead iterates
#' over the offending \emph{edges} and flags only the higher-scoring
#' endpoint per edge.  A genuine spike (with two violating edges) gets
#' flagged via either edge by virtue of its higher score and ends up in
#' the \code{is_outlier} set exactly once; its neighbours, with lower
#' scores, are spared.  Convergence is still guaranteed (the offending
#' edge is broken once one endpoint is removed) but may take more
#' iterations than the symmetric variant on dense error clusters.
#'
#' Tie-breaking: when \code{score_left == score_right}, the right
#' (later-in-time) endpoint is flagged.  This is arbitrary; users with
#' frequent ties should pre-jitter their scores.
#'
#' \code{aux_scores} are not modified or re-aligned during peel; the
#' user supplies a vector aligned to the original input rows and the
#' function indexes into it.  Score reliability is the user's
#' responsibility -- supplying noisy or anti-correlated scores can
#' produce worse results than the symmetric default.
#'
#' \strong{Cluster-outlier caveat.}  Asymmetric peel is designed for
#' \emph{1-fix spikes} where one endpoint of an offending edge is the
#' true outlier and the other is a clean neighbour.  On
#' \emph{coherent multi-fix clusters} (e.g. a long spoof segment that
#' is locally self-consistent and only violates \code{v_max} at the
#' jump-out / jump-back boundaries), the auxiliary score may not
#' discriminate between the spoof-interior boundary fix and the
#' clean-track boundary fix.  The peel may then walk inward from one
#' side only, or oscillate.  For known cluster-outlier datasets
#' (Argos PTT spoof, K02-style), the symmetric default is the
#' robust choice.
#'
#' @return The input \code{x} with added columns:
#'   \describe{
#'     \item{\code{is_outlier}}{Logical.  \code{TRUE} for fixes peeled
#'       at any iteration.}
#'     \item{\code{peel_iteration}}{Integer.  The iteration at which a
#'       fix was peeled, or \code{NA} if it survived.}
#'     \item{\code{step_speed}}{Numeric, m/s.  The outgoing step speed
#'       on the survivor sequence (\code{NA} where the fix was
#'       peeled or is at a track boundary).}
#'   }
#'   Attributes:
#'   \describe{
#'     \item{\code{v_max_used}}{The cap supplied.}
#'     \item{\code{n_peel_iterations}}{Integer, iterations performed.}
#'     \item{\code{converged}}{Logical.  \code{TRUE} if the peel
#'       finished before hitting \code{max_iter}.}
#'   }
#'
#' @seealso
#' \code{\link{mt_suggest_speed_cap}} for diagnostic inspection of
#' the speed distribution before choosing \code{v_max}.
#' \code{\link{mt_flag_speed_cap}} for one-shot flagging without
#' iteration.
#' \code{\link{mt_clean_track}} for the full per-fix + speed pipeline,
#' which uses \code{mt_peel_speed} internally when \code{v_max} is
#' supplied.
#'
#' @examples
#' \dontrun{
#' ## Eagle with a known spoof cluster; physiological cap ~30 m/s.
#' clean <- mt_peel_speed(eagle_track, v_max = 30)
#' summary(clean$is_outlier)
#' table(clean$peel_iteration, useNA = "ifany")
#' }
#'
#' @importFrom move2 mt_track_id mt_time
#' @importFrom sf st_coordinates st_is_longlat
#' @export
mt_peel_speed <- function(x, v_max, aux_scores = NULL, max_iter = 1000L,
                            remove = FALSE, silent = FALSE) {
  say <- function(...) if (!silent) message(...)

  ## ---- input validation ------------------------------------------
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  if (missing(v_max) || !is.numeric(v_max) || length(v_max) != 1L ||
      is.na(v_max) || v_max <= 0) {
    rlang::abort("`v_max` must be a positive scalar (m/s).",
                 class = "move2utils_mt_peel_speed_bad_v_max")
  }
  if (!is.null(aux_scores)) {
    if (!is.numeric(aux_scores) || length(aux_scores) != nrow(x)) {
      rlang::abort(
        sprintf("`aux_scores` must be a numeric vector of length nrow(x) = %d.",
                nrow(x)),
        class = "move2utils_mt_peel_speed_bad_aux_scores")
    }
  }
  max_iter <- as.integer(max_iter)
  if (is.na(max_iter) || max_iter < 1L) {
    rlang::abort("`max_iter` must be a positive integer.",
                 class = "move2utils_mt_peel_speed_bad_max_iter")
  }

  was_longlat <- isTRUE(sf::st_is_longlat(x))

  ## ---- common plumbing: multi-track dispatch + extraction + hygiene
  ## + .peel_fn_core + lift.  Peel stays in input CRS (Haversine on
  ## lon/lat, Euclidean on projected) for the same large-track reasons
  ## as detour. ----
  out <- .clean_track_dispatch(
    x,
    fn_core             = .peel_fn_core,
    fn_core_args        = list(was_longlat = was_longlat,
                                v_max       = v_max,
                                max_iter    = max_iter,
                                silent      = silent),
    per_track_args      = list(aux_scores  = aux_scores),
    need_time           = TRUE,
    hygiene_strict_time = FALSE,
    project_longlat     = FALSE,
    n_min               = 3L,
    n_min_severity      = "say",
    primitive_label     = "peel",
    silent              = silent)

  ## peel-specific attribute: v_max_used is just the input parameter
  attr(out, "v_max_used") <- v_max

  if (remove) out <- out[!out$is_outlier, ]
  out
}


## Raw-matrix entry point used by mt_clean_track to skip per-iter
## sf-class slicing.  Operates on (cc, t_s, active_idx, was_longlat,
## v_max, aux_scores, max_iter, silent).  Returns active-indexed
## result vectors (is_outlier, peel_iteration, step_speed) plus
## scalar attributes (n_peel_iterations, converged) lifted by the
## dispatcher.
##
## Stays in input CRS (Haversine on lon/lat, Euclidean on projected)
## to avoid AEQD projection cost on multi-million-fix tracks.
##
## @keywords internal
.peel_fn_core <- function(cc, t_s, active_idx, was_longlat,
                           v_max, aux_scores = NULL, max_iter = 1000L,
                           silent = TRUE) {

  ## Contract: active_idx sorted ascending with no duplicates.
  stopifnot(!is.unsorted(active_idx), !anyDuplicated(active_idx))

  say <- function(...) if (!silent) message(...)
  n   <- length(t_s)

  ## Initial keep mask: TRUE only at positions in active_idx.
  keep            <- rep(FALSE, n)
  keep[active_idx] <- TRUE

  peel_iter       <- rep(NA_integer_, n)
  final_speed_out <- rep(NA_real_, n)
  converged       <- FALSE
  it              <- 0L

  for (it in seq_len(max_iter)) {
    active <- which(keep)
    m <- length(active)
    if (m < 3L) {
      converged <- TRUE
      break
    }
    cc_k <- cc[active, , drop = FALSE]
    t_k  <- t_s[active]
    step_m <- if (was_longlat) {
      .haversine_pair(cc_k[-m, 1L], cc_k[-m, 2L],
                      cc_k[-1L, 1L], cc_k[-1L, 2L])
    } else {
      sqrt(diff(cc_k[, 1L])^2 + diff(cc_k[, 2L])^2)
    }
    dt         <- diff(t_k)
    step_speed <- ifelse(is.finite(dt) & dt > 0, step_m / dt, NA_real_)
    speed_in_per_fix  <- c(NA_real_, step_speed)
    speed_out_per_fix <- c(step_speed, NA_real_)

    if (is.null(aux_scores)) {
      ## Symmetric peel: any fix with max(in, out) > v_max is flagged.
      fix_speed <- pmax(speed_in_per_fix, speed_out_per_fix, na.rm = TRUE)
      violate <- is.finite(fix_speed) & fix_speed > v_max
      if (!any(violate)) {
        converged <- TRUE
        final_speed_out[active] <- speed_out_per_fix
        break
      }
      final_speed_out[active] <- speed_out_per_fix
      keep[active[violate]] <- FALSE
      peel_iter[active[violate]] <- it
    } else {
      ## Asymmetric peel: walk offending edges, flag the higher-
      ## scoring endpoint per edge.
      edge_violates <- is.finite(step_speed) & step_speed > v_max
      if (!any(edge_violates)) {
        converged <- TRUE
        final_speed_out[active] <- speed_out_per_fix
        break
      }
      e_idx       <- which(edge_violates)
      left_orig   <- active[e_idx]
      right_orig  <- active[e_idx + 1L]
      score_left  <- aux_scores[left_orig]
      score_right <- aux_scores[right_orig]
      flag_right_only <- !is.na(score_left) & !is.na(score_right) &
        score_right >= score_left
      flag_left_only  <- !is.na(score_left) & !is.na(score_right) &
        score_right < score_left
      flag_both <- is.na(score_left) | is.na(score_right)
      flagged <- unique(c(
        right_orig[flag_right_only],
        left_orig[flag_left_only],
        left_orig[flag_both],
        right_orig[flag_both]
      ))
      final_speed_out[active] <- speed_out_per_fix
      keep[flagged] <- FALSE
      peel_iter[flagged] <- it
    }
  }

  ## If max_iter exhausted, final speed pass on whatever survived
  if (!converged) {
    active <- which(keep)
    m <- length(active)
    if (m >= 2L) {
      cc_k <- cc[active, , drop = FALSE]
      t_k  <- t_s[active]
      step_m <- if (was_longlat) {
        .haversine_pair(cc_k[-m, 1L], cc_k[-m, 2L],
                        cc_k[-1L, 1L], cc_k[-1L, 2L])
      } else {
        sqrt(diff(cc_k[, 1L])^2 + diff(cc_k[, 2L])^2)
      }
      dt          <- diff(t_k)
      speed_out_k <- ifelse(is.finite(dt) & dt > 0, step_m / dt, NA_real_)
      final_speed_out[active] <- c(speed_out_k, NA_real_)
    }
  }

  ## NA out step_speed for peeled fixes (see original wrapper comment).
  is_out <- !keep & seq_len(n) %in% active_idx  # only flag within original active subset
  final_speed_out[is_out] <- NA_real_

  n_peel_iterations <- if (converged) max(0L, it - 1L) else max_iter

  say(sprintf(
    "Speed peel: v_max = %g m/s -- %d fix(es) peeled in %d iteration%s%s.",
    v_max, sum(is_out),
    n_peel_iterations,
    if (n_peel_iterations == 1L) "" else "s",
    if (converged) "" else " (hit max_iter)"))

  list(
    is_outlier        = is_out          [active_idx],
    peel_iteration    = peel_iter       [active_idx],
    step_speed        = final_speed_out [active_idx],
    n_peel_iterations = n_peel_iterations,
    converged         = converged
  )
}

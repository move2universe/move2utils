## Bridge-residual outlier detection for move2 objects.
##
## The movement-metric-based methods in mt_flag_outliers() detect
## outliers by extreme *transitions* in speed, angular velocity, and
## their auto-differences.  They are blind to outliers that sit in the
## middle of otherwise plausible movement metrics -- a location that is
## spatially impossible but whose speed-in and speed-out are compatible
## with the surrounding bouts.  On continuous-behaviour tracks (real
## animals as opposed to discretely-stated simulations) these block-
## shaped and drift-cluster errors are often dominant.
##
## This function closes that gap by treating each fix as a spatial
## observation to be checked against the Brownian bridge interpolated
## from its neighbours.  The bridge residual (deviation from the
## expected position given neighbours and the elapsed times) is
## normalised by the bridge width (sqrt(dt1 * dt2 / (dt1 + dt2))) to
## yield a gap-aware, scale-free score.  The same broken-stick and
## tail-decay inflection thresholding used by mt_flag_outliers() then
## identifies a natural break between outliers and the bulk.


#' Flag outliers via Brownian bridge residuals
#'
#' Detects GPS outliers by scoring each fix against the Brownian bridge
#' interpolated from its temporal neighbours.  A point sitting far from
#' the bridge mean relative to the bridge width is spatially implausible
#' given its neighbours and the local sampling cadence.  Complements
#' \code{\link{mt_flag_outliers}} by catching block-shaped errors and
#' faint near-cluster displacements that movement-metric methods miss.
#'
#' @details
#' For each fix \eqn{i} with temporal neighbours \eqn{i-1} and \eqn{i+1}:
#'
#' \enumerate{
#'   \item The Brownian-bridge mean at \eqn{i} is the time-weighted
#'         midpoint of the neighbours' coordinates.
#'   \item The bridge width is \eqn{\sqrt{\Delta t_1 \Delta t_2 /
#'         (\Delta t_1 + \Delta t_2)}} where \eqn{\Delta t_1,
#'         \Delta t_2} are the gaps to the previous and next fix.
#'   \item The score \eqn{\eta_i} is the Euclidean residual divided by
#'         the bridge width.  It is gap-aware by construction and does
#'         not depend on a local variance estimate (which would be
#'         inflated by the outlier itself -- the leverage problem).
#' }
#'
#' Outliers are identified by applying the broken-stick plus
#' tail-decay inflection threshold (same algorithm as
#' \code{\link{mt_flag_outliers}}) to \eqn{-\log \eta}.  Outliers sit
#' in the lower tail of \eqn{-\log \eta}, equivalently the upper tail
#' of \eqn{\eta}.
#'
#' Two post-processing steps clean the result:
#'
#' \itemize{
#'   \item \strong{Neighbour dedup.}  An outlier at fix \eqn{i}
#'         distorts the bridges centred on \eqn{i-1} and \eqn{i+1},
#'         inflating their \eqn{\eta} too.  For each run of
#'         consecutively-flagged indices, only the one with the
#'         largest \eqn{\eta} is kept.
#'   \item \strong{Iterative refinement.}  After flagging, residuals
#'         and thresholds are recomputed on the remaining points, and
#'         the process repeats until no new flags appear or
#'         \code{iterations} is reached.  This catches outliers that
#'         were masked by nearby ones on the first pass.
#' }
#'
#' The bridge-residual math is Euclidean and therefore requires metric
#' coordinates.  If the input is in longitude/latitude the function
#' auto-projects to a local azimuthal-equidistant CRS (centred on the
#' track's centroid) for the computation, and returns the result in
#' your original CRS.  \code{bridge_residual} and \code{bridge_width}
#' are always in metres regardless of the input CRS.
#'
#' @section Input preprocessing:
#' The bridge detector assumes each track is a time-sorted sequence
#' of finite-coordinate fixes with no duplicate timestamps. The
#' recommended preprocessing chain on a raw Movebank download is:
#'
#' \preformatted{
#' x <- mt_filter_gps_quality(x)              # sat/DOP/hacc + empty geoms
#' x <- move2::mt_filter_unique(x, "first")   # resolve duplicate times
#' x <- dplyr::arrange(x, mt_time(x))         # ensure time-sorted
#' x <- mt_flag_outliers_bridge(x)            # auto-projects if lon/lat
#' }
#'
#' Empty geometries and non-finite (\code{NA}) coordinates are a hard
#' error: there is no outlier to identify in a fix with no location, so
#' remove them upstream first (e.g. \code{x[!sf::st_is_empty(x), ]} or
#' \code{mt_filter_gps_quality(x)}). Rows with missing timestamps are
#' excluded from scoring (with an input-hygiene message) but retained in
#' the output. Duplicate or out-of-order timestamps are a hard error,
#' because their editorial resolution (which of two coincident fixes to
#' keep) is the user's decision, not this function's.
#'
#' @param x A \code{move2} object in any CRS (longitude/latitude is
#'   auto-projected to a local AEQD internally; the output is returned in
#'   the caller's CRS).  Single- or multi-track; multi-track inputs are
#'   processed per-individual.
#' @param method Character, one of \code{"combined"} (default),
#'   \code{"isotropic"}, or \code{"directional"}.  The residual
#'   decomposition gives two correlated but non-redundant scores per fix:
#'   \itemize{
#'     \item \eqn{\eta} — scalar residual / width (\code{"isotropic"});
#'           catches strong-omnidirectional outliers.
#'     \item \eqn{\eta_\perp} — orthogonal residual / width
#'           (\code{"directional"}); catches perpendicular-leverage
#'           outliers whose parallel component would otherwise dilute
#'           the scalar signal.  Concentrates on multipath errors and
#'           spoofing clusters, which drift across-track rather than
#'           along-track.
#'   }
#'   \code{"combined"} applies \code{threshold_type} independently to
#'   each score and flags a fix if \emph{either} score trips ---
#'   Pareto-dominates the single-score methods on synthetic
#'   ground-truth benchmarks and is the recommended default.
#'   \code{"isotropic"} and \code{"directional"} restrict flagging to a
#'   single score for interpretation or for reproducing the individual
#'   primitives in isolation.
#'
#'   Note: earlier versions named these options \code{"dBBMM"} and
#'   \code{"dBGB"} in reference to the dynamic Brownian-bridge
#'   framework.  The names were misleading --- this method uses the
#'   bridge-mean construction and the temporal width factor
#'   \eqn{\sqrt{\Delta t_1 \Delta t_2 / (\Delta t_1 + \Delta t_2)}}
#'   but deliberately does \emph{not} invoke the variance-estimation
#'   machinery with its window and margin that defines the real
#'   dBBMM / dBGB (see \code{\link{mt_dbbmm_variance}} and
#'   \code{\link{mt_dbgb_variance}} for those).  Omitting the variance
#'   is intentional: a locally-estimated variance is corrupted by the
#'   very outliers it would denominate (the leverage problem).
#' @param threshold_type Character, one of \code{"entropy"} (default)
#'   or \code{"gap"}.  \code{"entropy"} requires a real density valley
#'   between outliers and the bulk of \eqn{-\log \eta}; it returns no
#'   outliers on clean tracks.  \code{"gap"} uses the broken-stick and
#'   tail-decay inflection detector; it is more sensitive but can
#'   over-flag unimodal tails.
#' @param threshold Numeric.  Detection strictness.  Its meaning
#'   depends on \code{threshold_type}:
#'   \code{"entropy"} uses it as the maximum valley-to-peak density
#'   ratio (default 0.3); \code{"gap"} uses it as the break-size
#'   multiplier (default 3).  If \code{NULL}, a sensible default is
#'   chosen based on \code{threshold_type}.
#' @param location_error Per-fix observation-error prior, in metres
#'   (1-sigma horizontal).  Default \code{NULL} (disabled, current
#'   behaviour).  Accepts:
#'   \itemize{
#'     \item \code{NULL} -- no obs-error injection.
#'     \item a positive numeric scalar -- uniform sigma applied to
#'           every fix.  Useful when the device has a quoted nominal
#'           accuracy but no per-fix quality column.
#'     \item a numeric vector of length \code{nrow(x)} -- per-fix
#'           sigma already in metres.
#'     \item a single character string -- name of a column in \code{x}
#'           containing per-fix sigma in metres
#'           (e.g. \code{"eobs_horizontal_accuracy_estimate"}).
#'     \item the literal \code{"auto"} -- probes
#'           \code{eobs_horizontal_accuracy_estimate} first, then
#'           \code{argos_lc} (mapped via the standard CLS / Vincent
#'           et al. 2002 location-class table).  Falls back to no
#'           injection (with a message) when neither column is
#'           present.
#'   }
#'   When supplied, the bridge denominator at each fix is augmented
#'   with the variance contribution from its anchors only -- the
#'   target fix's own sigma deliberately does not enter, preserving
#'   leverage immunity.  The anchor variance is converted into
#'   bridge-width-equivalent units via an empirical residual-scale
#'   estimator \eqn{\hat S = \mathrm{median}(r_i^2 / w_i^2)} computed
#'   on the active mask in iteration 1.  See the diagnostic column
#'   \code{bridge_obs_inflation} for how strongly each fix's denominator
#'   was inflated.  The mechanism is most useful for Argos / mixed-mode
#'   tracks where per-fix accuracy varies by orders of magnitude;
#'   on uniform-quality GPS tracks it is typically a small correction.
#' @param residual_floor Numeric, non-negative.  Minimum absolute
#'   bridge residual (metres) for a rate-flagged fix to actually be
#'   flagged as an outlier.  Default \code{0} (disabled); the
#'   pre-2026 rate-only behaviour.  Set to a positive value (commonly
#'   5--25 m, or your device's nominal accuracy) to add a two-axis
#'   criterion: a fix is flagged only where both the rate score trips
#'   the threshold AND the absolute residual exceeds this floor.
#'   Recommended for real-world data with burst-mode sampling, where
#'   a few-metre GPS jitter over a 1-second dt produces a large eta
#'   that is not a physical outlier.  Not on by default so existing
#'   synthetic-ground-truth benchmarks (which can include sub-noise
#'   displacement outliers by construction) continue to pass.
#' @param iterations Integer.  Maximum number of refinement passes.
#'   Default 3 is HEURISTIC -- empirically the bulk of refinement
#'   gain happens in the first 2--3 iterations on the synthetic CPF
#'   benchmark.  Iteration stops early when no new points are flagged.
#'   Plausible range: 2--5.  Higher values increase runtime without
#'   meaningful F1 gain on the benchmark cohort.
#' @param dedup_neighbours Logical.  If \code{TRUE} (default), remove
#'   neighbour-smearing by keeping only the peak \eqn{\eta} in each
#'   consecutive run of flagged indices.
#' @param pool_by Optional character vector of length 1 or 2 naming
#'   column(s) in \code{mt_track_data(x)}.  Length 1: single column
#'   used as both fit set and operating unit.  Length 2:
#'   \code{c(outer, inner)} where \code{outer} names the fit-source
#'   column (the union of its events supplies the entropy / gap
#'   break distribution on \code{bridge_eta} and, for
#'   \code{"directional"}/\code{"combined"}, \code{bridge_eta_perp})
#'   and \code{inner} names the operating unit (within which pool-
#'   added flags are unioned).  Length 2 requires strict nesting:
#'   every distinct \code{inner} value must map to exactly one
#'   \code{outer} value.  Length \eqn{> 2} is rejected (pool_by has
#'   exactly two semantic roles).  Pool flags union into
#'   \code{is_outlier} -- additive, never un-flags what per-track
#'   iteration caught.  The \code{residual_floor} gate is respected
#'   at the pool level (matches the per-track gate).  No dedup at
#'   the pool level (per-track dedup already ran inside the
#'   iteration).  \code{NULL} (default) preserves per-track
#'   behaviour byte-identically.  See \code{?mt_clean_track} for
#'   the orchestrator-level walkthrough.
#' @param plot Logical.  If \code{TRUE} (default), produce a
#'   diagnostic plot (sorted \eqn{\log \eta} with break + map of flags).
#' @param remove Logical.  If \code{TRUE}, return the object with
#'   flagged rows removed.  Default \code{FALSE}.
#' @param silent Logical.  If \code{FALSE} (default) the function
#'   prints a brief running narration: per-iteration break and flag
#'   counts, input-hygiene notes, projection messages, and a final
#'   summary line.  Set \code{TRUE} to suppress.  Errors and
#'   warnings are always shown.
#'
#' @return A \code{move2} object with columns added:
#'   \describe{
#'     \item{\code{bridge_residual}}{Euclidean deviation from the
#'           bridge mean, in metres.}
#'     \item{\code{bridge_width}}{The bridge-width factor
#'           \eqn{\sqrt{\Delta t_1 \Delta t_2 / (\Delta t_1 +
#'           \Delta t_2)}}, in \eqn{\sqrt{s}}.}
#'     \item{\code{bridge_eta}}{The normalised bridge score used for
#'           flagging.  For \code{method = "isotropic"} this is
#'           \eqn{\eta_i = r_i / w_i}.  For
#'           \code{method = "directional"} it is the orthogonal component
#'           \eqn{\eta_{\perp,i} = r_{\perp,i} / w_i}.}
#'     \item{\code{bridge_eta_para}, \code{bridge_eta_perp}}{Gap-
#'           normalised parallel and
#'           orthogonal residuals.  Useful for inspecting whether a
#'           flagged point is primarily an along-track or
#'           perpendicular anomaly.}
#'     \item{\code{bridge_obs_inflation}}{Ratio of effective to
#'           geometric bridge width \eqn{w_{\text{eff}}/w}.  Equal to
#'           \code{1} when \code{location_error} is \code{NULL} or all
#'           anchors have unknown sigma; \eqn{>1} where anchor
#'           obs-error contributed meaningfully.}
#'     \item{\code{bridge_percentile}}{Empirical percentile of
#'           \eqn{-\log \eta} (0 = most extreme outlier).}
#'     \item{\code{bridge_iteration}}{Integer iteration at which the
#'           point was flagged, or \code{NA} if unflagged.}
#'     \item{\code{loglr_bridge}}{Signed log-likelihood-ratio of outlier
#'       vs not, in nat units; \code{> 0} where flagged, \code{< 0} below
#'       the boundary, \code{NA} where the detector abstains.  Same
#'       convention as \code{loglr_prob}; see
#'       \code{DESIGN_evidence_accumulation.md}.}
#'     \item{\code{is_outlier}}{Logical flag.}
#'     \item{\code{flagged_by_bridge}}{Logical. Same as
#'       \code{is_outlier} for this primitive; named for parity with the
#'       other primitives' output schema (so a track scored by several
#'       primitives can be voted by \code{\link{mt_flag_consensus}}).}
#'     \item{\code{is_na_prob}}{\code{TRUE} where \eqn{\eta} is
#'           undefined (track endpoints, missing neighbours).}
#'   }
#'   If \code{remove = TRUE}, the flagged rows are dropped.
#'
#' @examples
#' \dontrun{
#' library(move2)
#' syn <- mt_read(system.file("extdata", "synthetic_tracks.csv.gz",
#'                             package = "move2utils"))
#' ## any CRS works -- lon/lat is auto-projected internally
#' result <- mt_flag_outliers_bridge(syn)
#' }
#'
#' @seealso \code{\link{mt_flag_outliers}} for movement-metric
#'   detection; \code{\link{mt_combined_outliers}} for composing
#'   detectors via majority vote.
#'
#' @references
#' Safi, K. (2026). Self-thresholding hierarchical
#' outlier-detection for animal movement tracks. Companion paper to
#' the \pkg{move2utils} R package. bioRxiv preprint, submitted to
#' Methods in Ecology and Evolution. \doi{10.64898/2026.07.11.737894}
#'
#' @importFrom move2 mt_time mt_track_id mt_n_tracks
#' @importFrom sf st_coordinates st_is_longlat
#' @importFrom graphics par plot points lines abline legend
#' @importFrom grDevices adjustcolor
#' @export
mt_flag_outliers_bridge <- function(x,
                                     method = c("combined", "isotropic", "directional"),
                                     threshold_type = c("entropy", "gap"),
                                     threshold = NULL,
                                     location_error = NULL,
                                     residual_floor = 0,
                                     iterations = 3,
                                     dedup_neighbours = TRUE,
                                     pool_by = NULL,
                                     plot = TRUE,
                                     remove = FALSE,
                                     silent = FALSE) {
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  method         <- match.arg(method)
  ## brief narrator helper -- suppressed under silent = TRUE
  say <- .say(silent)
  threshold_type <- match.arg(threshold_type)
  if (!is.null(pool_by)) {
    ## `.resolve_pool_groups` validates shape, columns, and (for
    ## length-2 input) the strict-nesting requirement.
    invisible(.resolve_pool_groups(x, pool_by, silent = TRUE))
  }
  ## Default thresholds are HEURISTIC and live on the leaf detectors
  ## as their formal defaults:
  ##
  ##   entropy = 0.3 (.entropy_threshold_lower$threshold): Maximum
  ##     valley-to-peak density ratio in the log-residual KDE. Unified
  ##     package-wide entropy default (2026-05-09); validated by the
  ##     2026-05-06 Raven sensitivity sweep as the only level inside
  ##     the strict stability window for cohort flag rate. Cascade
  ##     flag rate is statistically insensitive in [0.3, 0.7] (K-W
  ##     p = 0.076). Plausible range: 0.3--0.7.
  ##
  ##   gap = 3 (.gap_threshold_lower$threshold): Break-size multiplier
  ##     for the broken-stick detector. The "3-sigma" convention
  ##     applied to log-residual gaps. Plausible range: 2--5.
  ##
  ## `threshold = NULL` (the default here) defers to the leaf formal,
  ## the single source of truth.  A non-NULL scalar overrides the leaf.
  ## See audits/2026-05-25-parameter-propagation/findings.md §1.
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        is.na(threshold) || threshold <= 0) {
      rlang::abort("`threshold` must be a positive scalar or NULL.",
                   class = "move2utils_mt_flag_outliers_bridge_bad_threshold")
    }
  }
  if (!is.numeric(residual_floor) || length(residual_floor) != 1 ||
      is.na(residual_floor) || residual_floor < 0) {
    rlang::abort("`residual_floor` must be a non-negative scalar (metres).",
                 class = "move2utils_mt_flag_outliers_bridge_bad_residual_floor")
  }

  ## ---- resolve obs-error specification ---------------------------
  ## Resolved on the original `x` so column lookups (e.g. `"auto"` or a
  ## column name) work regardless of CRS / sliced state.  The result
  ## travels through `.clean_track_dispatch`'s per-track slicing via
  ## `per_track_args`.
  obs_sigma <- .resolve_location_error(location_error, x, nrow(x))
  if (!is.null(obs_sigma) && all(is.na(obs_sigma))) {
    say("  location_error: all per-fix sigmas are NA; injection disabled.")
    obs_sigma <- NULL
  }

  ## ---- common plumbing: multi-track dispatch + projection + extraction
  ## + hygiene + .bridge_fn_core + full-length lift + CRS restore ----
  out <- .clean_track_dispatch(
    x,
    fn_core             = .bridge_fn_core,
    fn_core_args        = list(method           = method,
                                threshold_type   = threshold_type,
                                threshold        = threshold,
                                residual_floor   = residual_floor,
                                iterations       = iterations,
                                dedup_neighbours = dedup_neighbours,
                                silent           = silent),
    per_track_args      = list(location_error   = obs_sigma),
    need_time           = TRUE,
    hygiene_strict_time = TRUE,
    n_min               = 10L,
    primitive_label     = "bridge",
    silent              = silent)

  ## ---- bridge-specific post-processing -----------------------------
  ## `is_na_prob` lives downstream of the lift: non-active rows have
  ## `bridge_eta = NA` (no probability computed), and active rows with
  ## `bridge_eta` non-finite or <= 0 also count as missing-probability.
  out$is_na_prob <- !is.finite(out$bridge_eta) | out$bridge_eta <= 0
  attr(out, "bridge_method") <- method

  ## ---- pool_by union: refit threshold(s) per group and union flags ----
  ## Approach (ii): per-track dispatch runs with track-local iteration
  ## (preserving today's behaviour); then a pool threshold is fit per
  ## group from the union of converged bridge_eta / bridge_eta_perp
  ## values and the pool flags are unioned into is_outlier.  Pool
  ## flags are ADDITIVE only -- never un-flag what the per-track
  ## iteration caught.  Approximation note: per-track iteration uses
  ## track-local thresholds inside its inner loop; pool threshold is
  ## applied at the end to the final eta vector.  For the one-pass
  ## primitives (detour, speed_cap) this is exact; for bridge it is
  ## an approximation that converges to the exact pooled iteration
  ## when track-local and pool thresholds agree (which they do when
  ## per-track distribution is representative of the pooled one).
  if (!is.null(pool_by)) {
    out <- .bridge_pool_union(out, x, pool_by,
                                method         = method,
                                threshold_type = threshold_type,
                                threshold      = threshold,
                                residual_floor = residual_floor,
                                silent         = silent)
  }

  ## log-LR emission (evidence currency; additive, is_outlier untouched)
  out <- .attach_bridge_loglr(out, method)

  if (plot)   .plot_bridge_outliers(out)
  if (remove) out <- out[!out$is_outlier, ]
  out
}

## Bridge log-LR.  `eta` (= residual / width) is a standardized residual
## magnitude |Z|, higher = more outlier, so the Gaussian surprisal is
## eta^2 / 2 nats -- NOT -log(eta).  The relevant axis is eta (isotropic),
## eta_perp (directional), or the larger of the two (combined flags on
## either break -> the more surprising axis).  Competence (bridge_width vs
## local use-scale) is layered on in step 3b.
## @keywords internal
.attach_bridge_loglr <- function(obj, method = "combined") {
  ## flag-source column: makes the standalone primitive votable by
  ## mt_flag_consensus (identical to is_outlier for a single detector;
  ## travels with loglr_bridge through the same back-attach path).
  obj$flagged_by_bridge <- obj$is_outlier
  eta <- obj$bridge_eta
  if (is.null(eta)) return(obj)
  z <- switch(method,
    isotropic   = eta,
    directional = obj$bridge_eta_perp,
    combined    = pmax(eta, obj$bridge_eta_perp, na.rm = TRUE),
    eta)
  neglogp <- 0.5 * z^2
  if (!is.null(obj$is_na_prob)) {
    neglogp[as.logical(obj$is_na_prob) %in% TRUE] <- NA_real_
  }
  obj$loglr_bridge <- .loglr_grouped(neglogp, obj$is_outlier,
                                     move2::mt_track_id(obj))
  obj
}


## Pool-union closure for mt_flag_outliers_bridge.
##
## For each pool group, fits one (or two, for method="combined")
## thresholds on the union of converged bridge_eta / bridge_eta_perp
## values across the group's tracks.  Applies the pool break per
## track, respecting the same residual_floor gate the per-track pass
## uses.  Unions into is_outlier.
##
## Method dispatch:
##   isotropic   -> pool break on neg_log(bridge_eta)        only
##   directional -> pool break on neg_log(bridge_eta_perp)   only
##   combined    -> both pool breaks; flags union pre-dedup
##
## Dedup: no dedup at the pool level.  Per-track dedup already ran
## inside .bridge_fn_core's iteration; pool flags are additive on
## top of dedup-applied converged per-track output.
##
## @keywords internal
.bridge_pool_union <- function(out, x, pool_by, method,
                                threshold_type, threshold,
                                residual_floor, silent = FALSE) {

  ## NULL threshold defers to the leaf formal (single source of truth).
  ## Local threshold-fit helper matching .bridge_fn_core's thr_fn.
  thr_fn <- function(v) {
    leaf <- switch(threshold_type,
                    entropy = .entropy_threshold_lower,
                    gap     = .gap_threshold_lower)
    if (is.null(threshold)) leaf(v)
    else                    leaf(v, threshold = threshold)
  }

  neg_log <- function(v) {
    ok  <- is.finite(v) & v > 0
    nl  <- rep(NA_real_, length(v))
    nl[ok] <- -log(v[ok])
    nl
  }

  pool_step <- function(fit_idx, apply_idx) {
    ## Fit threshold(s) on the OUTER group; apply to the INNER group.
    eta_fit   <- out$bridge_eta[fit_idx]
    eta_p_fit <- out$bridge_eta_perp[fit_idx]
    eta_g     <- out$bridge_eta[apply_idx]
    eta_p_g   <- out$bridge_eta_perp[apply_idx]
    resid_g   <- out$bridge_residual[apply_idx]
    r_above   <- is.finite(resid_g) & resid_g > residual_floor

    ## Threshold-fit utilities return a break value implicitly tied
    ## to the input distribution.  We need the break value, then
    ## evaluate it against the inner-group eta values explicitly.
    apply_with <- function(eta_apply, fitted) {
      ## Mirror `.entropy_threshold_lower` / `.gap_threshold_lower`:
      ## a value is flagged when its neg-log score is strictly below
      ## the fitted break.
      bv <- fitted$break_value
      if (!is.numeric(bv) || length(bv) != 1L || is.na(bv)) {
        return(rep(FALSE, length(eta_apply)))
      }
      nl <- neg_log(eta_apply)
      !is.na(nl) & nl < bv
    }

    if (method == "isotropic") {
      thr_e <- thr_fn(neg_log(eta_fit))
      flagged <- apply_with(eta_g, thr_e) & r_above
    } else if (method == "directional") {
      thr_p <- thr_fn(neg_log(eta_p_fit))
      flagged <- apply_with(eta_p_g, thr_p) & r_above
    } else {  # combined
      thr_e <- thr_fn(neg_log(eta_fit))
      thr_p <- thr_fn(neg_log(eta_p_fit))
      flagged <- (apply_with(eta_g,   thr_e) |
                  apply_with(eta_p_g, thr_p)) & r_above
    }
    flagged[is.na(flagged)] <- FALSE
    if (!silent) {
      n_new <- sum(flagged & !out$is_outlier[apply_idx], na.rm = TRUE)
      message(sprintf(
        "  pool_by[%s]: %d new fix(es) flagged across %d event(s) (method = %s).",
        paste(pool_by, collapse = ","), n_new, length(apply_idx), method))
    }
    flagged
  }

  .apply_pool_union(out, x, pool_by, pool_step,
                    flag_cols = "is_outlier",
                    silent = silent)
}


## ---- internal entry point used by mt_clean_track -----------------------

#' Bridge-residual scoring, raw-matrix entry point.
#'
#' Internal counterpart to [mt_flag_outliers_bridge()].  Operates on
#' raw coordinate / time matrices rather than a `move2` object so that
#' [mt_clean_track()] can cache `(cc, t_s)` once and call the math
#' directly across iterations without re-extracting them per call.
#'
#' @param cc Numeric matrix, two columns, full-track length.  Caller is
#'   responsible for projection (Euclidean math expects a metric CRS).
#' @param t_s Numeric vector, full-track length.  Time in seconds.
#' @param active_idx Integer vector.  Row indices into `cc` and `t_s`
#'   that participate in this call.  The caller is responsible for
#'   excluding rows with non-finite coordinates / times and for
#'   ensuring no duplicate or out-of-order timestamps among the
#'   active subset.
#' @param method,threshold_type,threshold,residual_floor,iterations,dedup_neighbours
#'   As in [mt_flag_outliers_bridge()].
#' @param location_error Either `NULL` or a numeric vector of per-fix
#'   horizontal 1-sigma values (metres), **full-track length** (the
#'   helper indexes into it via `active_idx`).  `NA` anchors contribute
#'   zero obs-error variance.
#' @param silent Logical; suppress per-iteration narration.
#'
#' @return A list of vectors aligned to `active_idx`:
#' \describe{
#'   \item{`is_outlier`}{logical}
#'   \item{`bridge_residual`,`bridge_width`,`bridge_eta`,`bridge_eta_para`,`bridge_eta_perp`,`bridge_obs_inflation`,`bridge_percentile`}{numeric}
#'   \item{`bridge_iteration`}{integer}
#'   \item{`S_hat`}{scalar (`NULL` if obs-error injection disabled)}
#' }
#'
#' Caller lifts these to full-length output columns.
#'
#' @keywords internal
.bridge_fn_core <- function(cc, t_s, active_idx,
                            method            = "combined",
                            threshold_type    = "entropy",
                            threshold         = NULL,
                            location_error    = NULL,
                            residual_floor    = 0,
                            iterations        = 3L,
                            dedup_neighbours  = TRUE,
                            silent            = TRUE) {

  ## Contract: active_idx must be sorted ascending with no duplicates.
  ## The cascade + wrappers guarantee this via which(), but the math
  ## (diff(t_s[active_idx]), bridge weighting, dedup) silently corrupts
  ## under a shuffled or duplicated active_idx -- enforce explicitly.
  stopifnot(!is.unsorted(active_idx), !anyDuplicated(active_idx))

  k   <- length(active_idx)
  say <- function(...) if (!silent) message(...)

  ## Length-k internal storage.  This preserves the historical "dedup
  ## sees the active subset" semantics from the era when the cascade
  ## sliced its input (`active_x = x[active_idx, ]`) before calling the
  ## bridge wrapper -- already-flagged positions sit OUTSIDE the array,
  ## so two new flags separated by a previously-flagged fix are
  ## consecutive in active-subset space and `.dedup_consecutive`
  ## eliminates one (the correct cascade behaviour).  Length-n
  ## storage would have left previously-flagged positions in the array
  ## as FALSE, breaking that semantic and over-flagging on
  ## cascade-active-only inputs.
  active_local <- rep(TRUE, k)

  residual_out   <- rep(NA_real_,    k)
  width_out      <- rep(NA_real_,    k)
  eta_out        <- rep(NA_real_,    k)
  eta_para_out   <- rep(NA_real_,    k)
  eta_perp_out   <- rep(NA_real_,    k)
  inflation_out  <- rep(NA_real_,    k)
  pct_out        <- rep(NA_real_,    k)
  iter_out       <- rep(NA_integer_, k)
  is_outlier     <- rep(FALSE,       k)

  obs_sigma <- location_error  # full-length-n numeric or NULL
  S_hat     <- NULL

  say(sprintf("Running bridge-residual detection (method = %s) on %d locations...",
              method, k))

  for (iter in seq_len(iterations)) {
    local_alive <- which(active_local)
    if (length(local_alive) < 10) break

    ## Global indices into cc/t_s for the math helper.
    global_alive <- active_idx[local_alive]

    if (iter == 1L && !is.null(obs_sigma) && is.null(S_hat)) {
      br0 <- .compute_bridge_residuals_dbgb(cc, t_s, global_alive)
      rw2 <- (br0$residual / br0$width)^2
      S_hat <- stats::median(rw2, na.rm = TRUE)
      if (!is.finite(S_hat) || S_hat <= 0) {
        rlang::warn(
          "Could not estimate residual scale; location_error injection disabled.",
          class = "move2utils_mt_flag_outliers_bridge_residual_scale_failed")
        obs_sigma <- NULL
        S_hat     <- NULL
      } else {
        say(sprintf(
          "  location_error: residual scale S_hat = %.4g m^2/s; injecting anchor obs-error.",
          S_hat))
      }
    }

    br <- .compute_bridge_residuals_dbgb(cc, t_s, global_alive,
                                          sigma = obs_sigma, S_hat = S_hat)

    ## Write into local (k-sized) storage at the still-alive positions.
    residual_out [local_alive] <- br$residual
    width_out    [local_alive] <- br$width
    eta_out      [local_alive] <- br$eta
    eta_para_out [local_alive] <- br$eta_para
    eta_perp_out [local_alive] <- br$eta_perp
    inflation_out[local_alive] <- br$inflation

    neg_log <- function(v) {
      ok <- is.finite(v) & v > 0
      out <- rep(NA_real_, length(v))
      out[ok] <- -log(v[ok])
      out
    }
    thr_fn <- function(x) {
      leaf <- switch(threshold_type,
                      entropy = .entropy_threshold_lower,
                      gap     = .gap_threshold_lower)
      if (is.null(threshold)) leaf(x)
      else                    leaf(x, threshold = threshold)
    }

    thr_eta  <- thr_fn(neg_log(br$eta))
    thr_perp <- thr_fn(neg_log(br$eta_perp))

    pct_primary <- switch(method,
      isotropic   = thr_eta$percentile,
      directional = thr_perp$percentile,
      combined    = thr_eta$percentile)
    pct_out[local_alive] <- pct_primary

    r_above_floor <- is.finite(br$residual) & br$residual > residual_floor

    flagged_local_alive <- switch(method,
      isotropic   = thr_eta$is_outlier  & r_above_floor,
      directional = thr_perp$is_outlier & r_above_floor,
      combined    = (thr_eta$is_outlier | thr_perp$is_outlier) & r_above_floor)

    break_value <- switch(method,
      isotropic   = thr_eta$break_value,
      directional = thr_perp$break_value,
      combined    = NA_real_)

    if (method != "combined" && is.na(break_value)) {
      say(sprintf("  Iter %d: no break found (%s); stopping.",
                  iter, threshold_type))
      break
    }
    if (method == "combined" &&
        is.na(thr_eta$break_value) && is.na(thr_perp$break_value)) {
      say(sprintf("  Iter %d: no break on either score (%s); stopping.",
                  iter, threshold_type))
      break
    }

    n_flagged <- sum(flagged_local_alive, na.rm = TRUE)
    if (n_flagged == 0) {
      say(sprintf("  Iter %d: no outliers.", iter))
      break
    }

    ## Lift the alive-mask flags to local (k-sized) and dedup there.
    ## Dedup over length-k = "consecutive among currently-alive fixes",
    ## which is exactly the OLD wrapper's slice semantics when the
    ## cascade passed `active_x = x[active_idx, ]`.
    if (method == "combined") {
      gE <- rep(FALSE, k); gE[local_alive] <- thr_eta$is_outlier  & r_above_floor
      gP <- rep(FALSE, k); gP[local_alive] <- thr_perp$is_outlier & r_above_floor
      if (dedup_neighbours) {
        if (sum(gE) > 1) gE <- .dedup_consecutive(gE, eta_out)
        if (sum(gP) > 1) gP <- .dedup_consecutive(gP, eta_perp_out)
      }
      flagged_local_full <- gE | gP
    } else {
      flagged_local_full <- rep(FALSE, k)
      flagged_local_full[local_alive] <- flagged_local_alive
      dedup_score <- switch(method,
        isotropic   = eta_out,
        directional = eta_perp_out)
      if (dedup_neighbours && n_flagged > 1) {
        flagged_local_full <- .dedup_consecutive(flagged_local_full,
                                                  dedup_score)
      }
    }

    n_added <- sum(flagged_local_full & !is_outlier)
    if (n_added == 0) {
      say(sprintf("  Iter %d: %d flagged but all already-known; stopping.",
                  iter, n_flagged))
      break
    }

    iter_out[flagged_local_full & !is_outlier] <- iter
    is_outlier  [flagged_local_full] <- TRUE
    active_local[flagged_local_full] <- FALSE

    msg_break <- switch(method,
      isotropic   = sprintf("break at eta = %.2f",
                            exp(-thr_eta$break_value)),
      directional = sprintf("break at eta_perp = %.2f",
                            exp(-thr_perp$break_value)),
      combined    = sprintf("break at eta = %s / eta_perp = %s",
        if (is.na(thr_eta$break_value))  "-" else sprintf("%.2f", exp(-thr_eta$break_value)),
        if (is.na(thr_perp$break_value)) "-" else sprintf("%.2f", exp(-thr_perp$break_value))))
    say(sprintf("  Iter %d: flagged %d (%s).",
                iter, n_added, msg_break))
  }

  list(
    is_outlier           = is_outlier,      # already length-k
    bridge_residual      = residual_out,
    bridge_width         = width_out,
    bridge_eta           = eta_out,
    bridge_eta_para      = eta_para_out,
    bridge_eta_perp      = eta_perp_out,
    bridge_obs_inflation = inflation_out,
    bridge_percentile    = pct_out,
    bridge_iteration     = iter_out,
    bridge_S_hat         = S_hat
  )
}


## ---- helpers -----------------------------------------------------------


## Compute the directional (parallel / orthogonal) bridge residual for
## the active subsequence of a track.  Decomposition is relative to the
## local travel axis, defined as the direction from the bridge mean to
## the next active fix.  The orthogonal residual is the part of the
## displacement perpendicular to that axis -- the signature of
## multipath errors and spoofing clusters, which tend to drift
## sideways relative to true movement.
##
## When `sigma` (per-fix horizontal 1-sigma in m) and `S_hat` (the
## empirical residual scale, m^2/s; see mt_flag_outliers_bridge for
## construction) are both supplied, the bridge denominator is
## inflated to absorb anchor observation-error variance:
##
##   w_eff^2 = w^2 + (w_prev^2 sigma_prev^2 + w_next^2 sigma_next^2) / S_hat
##
## where w_prev = dt2/(dt1+dt2) and w_next = dt1/(dt1+dt2) match the
## bridge-mean weights.  Only anchor variances enter -- the target
## fix's sigma deliberately does not, preserving leverage immunity.
##
## Returns a list:
##   residual  : total Euclidean residual (m), same as the isotropic version
##   width     : sqrt(dt1 * dt2 / (dt1 + dt2))  -- geometric, sqrt(s)
##   width_eff : effective width incorporating anchor obs-error.
##               Equal to `width` when sigma/S_hat not supplied.
##   inflation : width_eff / width (>= 1; 1 when sigma/S_hat absent
##               or all anchors have sigma = 0)
##   eta_para  : parallel component magnitude / width_eff
##   eta_perp  : orthogonal component magnitude / width_eff
##   eta       : scalar (isotropic) residual / width_eff
##
## Scale normalisation: pure geometric in the no-obs-error case (same
## trick as the original isotropic version).  We divide by
## bridge_width rather than by a directional variance estimate from
## mt_dbgb_variance(), because that variance would be inflated by the
## very outliers we are trying to detect (the leverage problem).
## That is why this outlier method, despite being inspired by the
## dBBMM / dBGB bridge constructions, does NOT invoke their
## variance-estimation machinery.  The optional obs-error injection
## introduces a per-fix prior at the anchors only, which preserves
## that property.
##
## @keywords internal
.compute_bridge_residuals_dbgb <- function(cc, t, active_idx,
                                            sigma = NULL, S_hat = NULL) {
  nk <- length(active_idx)
  residual  <- rep(NA_real_, nk)
  width     <- rep(NA_real_, nk)
  width_eff <- rep(NA_real_, nk)
  inflation <- rep(NA_real_, nk)
  eta_para  <- rep(NA_real_, nk)
  eta_perp  <- rep(NA_real_, nk)

  if (nk < 3) {
    return(list(residual  = residual, width = width,
                width_eff = width_eff, inflation = inflation,
                eta_para  = eta_para, eta_perp  = eta_perp,
                eta = eta_perp))
  }

  cc_a <- cc[active_idx, , drop = FALSE]
  t_a  <- t[active_idx]

  dt1 <- diff(t_a)[-(nk - 1L)]
  dt2 <- diff(t_a)[-1L]

  ## bridge-mean weights: w_prev = dt2/(dt1+dt2), w_next = dt1/(dt1+dt2)
  w_next <- dt1 / (dt1 + dt2)
  w_prev <- 1 - w_next

  mu_x <- w_prev * cc_a[seq_len(nk - 2L), 1] + w_next * cc_a[3:nk, 1]
  mu_y <- w_prev * cc_a[seq_len(nk - 2L), 2] + w_next * cc_a[3:nk, 2]

  ## decomposition axis: unit vector from mu towards the NEXT active
  ## fix (cc_a[3:nk]).  This is the local travel direction.  Using
  ## mu->next (not prev->next) keeps the reference consistent with the
  ## .delta_para_orth convention in dbgb_variance.R.
  ax_x <- cc_a[3:nk, 1] - mu_x
  ax_y <- cc_a[3:nk, 2] - mu_y
  axn  <- sqrt(ax_x * ax_x + ax_y * ax_y)

  ## residual vector (observed - mu) at fix i
  rx <- cc_a[2:(nk - 1L), 1] - mu_x
  ry <- cc_a[2:(nk - 1L), 2] - mu_y
  rn <- sqrt(rx * rx + ry * ry)

  ## signed projection onto travel axis. is.finite() guards against
  ## NaN axn values (which would yield NA in has_axis and blow up the
  ## subscript assignment below). NaN propagates here when an active
  ## neighbour has a non-finite coord or when a triple of identical
  ## timestamps makes the time-weighting degenerate; those rows get
  ## NA eta via the `bad` filter further down.
  has_axis <- is.finite(axn) & axn > 0
  rp <- rep(NA_real_, nk - 2L)   # parallel magnitude
  ro <- rep(NA_real_, nk - 2L)   # orthogonal magnitude
  rp[has_axis] <-
    abs(rx[has_axis] * ax_x[has_axis] + ry[has_axis] * ax_y[has_axis]) /
    axn[has_axis]
  ## orthogonal magnitude: |residual_vec x axis_unit|, or via Pythagoras
  ro[has_axis] <- sqrt(pmax(rn[has_axis] ^ 2 - rp[has_axis] ^ 2, 0))
  ## when no axis direction (two coincident neighbours), fall back to
  ## isotropic: split total equally between para/orth via sqrt(2)
  rp[!has_axis] <- rn[!has_axis] / sqrt(2)
  ro[!has_axis] <- rn[!has_axis] / sqrt(2)

  bad <- !is.finite(dt1) | !is.finite(dt2) | dt1 <= 0 | dt2 <= 0
  rn[bad] <- NA_real_
  rp[bad] <- NA_real_
  ro[bad] <- NA_real_

  bw  <- sqrt(dt1 * dt2 / (dt1 + dt2))
  bw[bad] <- NA_real_

  ## effective bridge width (units of sqrt(s)): geometric in the
  ## no-obs-error case, anchor-error-augmented otherwise.  S_hat is
  ## the residual scale estimated by the caller from a first
  ## uncorrected pass (median of (residual / width)^2).
  use_obs <- !is.null(sigma) && !is.null(S_hat) &&
             is.finite(S_hat) && S_hat > 0
  if (use_obs) {
    sig_a    <- sigma[active_idx]
    sig_prev <- sig_a[seq_len(nk - 2L)]
    sig_next <- sig_a[3:nk]
    sig_prev[!is.finite(sig_prev)] <- 0   # NA anchor -> no contribution
    sig_next[!is.finite(sig_next)] <- 0
    obs_term <- (w_prev ^ 2) * (sig_prev ^ 2) +
                (w_next ^ 2) * (sig_next ^ 2)         # m^2
    bw2_eff  <- (dt1 * dt2 / (dt1 + dt2)) + obs_term / S_hat
    bw_eff   <- sqrt(bw2_eff)
    bw_eff[bad] <- NA_real_
  } else {
    bw_eff <- bw
  }

  residual[2:(nk - 1L)]  <- rn
  width[2:(nk - 1L)]     <- bw
  width_eff[2:(nk - 1L)] <- bw_eff
  inflation[2:(nk - 1L)] <- bw_eff / bw
  eta_para[2:(nk - 1L)]  <- rp / bw_eff
  eta_perp[2:(nk - 1L)]  <- ro / bw_eff

  ## scalar eta = ||r|| / w_eff  (isotropic; equals r/w when no obs-error)
  eta <- residual / width_eff

  list(residual  = residual,
       width     = width,
       width_eff = width_eff,
       inflation = inflation,
       eta       = eta,        # isotropic (scalar magnitude)
       eta_para  = eta_para,   # parallel (diagnostic only)
       eta_perp  = eta_perp)   # directional (orthogonal component)
}


## For each run of consecutively-flagged indices, keep only the one
## with the largest eta.  Suppresses neighbour-smearing: an outlier at
## i inflates eta at i-1 and i+1 through the bridge mean.
##
## @keywords internal
.dedup_consecutive <- function(flagged, eta) {
  idx <- which(flagged)
  if (length(idx) < 2) return(flagged)

  ## Group into runs of consecutive indices, then for each run keep
  ## the index with the largest eta. Vectorised via tapply over the
  ## run_id -- avoids a per-run for-loop that scaled poorly when a
  ## track had many small runs.
  run_id   <- cumsum(c(1L, diff(idx) > 1L))
  eta_safe <- ifelse(is.finite(eta[idx]), eta[idx], -Inf)
  keepers  <- unname(tapply(seq_along(idx), run_id,
                             function(i) idx[i[which.max(eta_safe[i])]]))
  out <- rep(FALSE, length(flagged))
  out[as.integer(keepers)] <- TRUE
  out
}


## Diagnostic plot for bridge outlier results.
##
## Left panel: directional decomposition of the residual on the
## (log eta_para, log eta_perp) plane.  Diagonal y = x separates
## isotropic jitter (on the diagonal) from along-track errors
## (below) and across-track errors (above).  Right panel: map.
##
## @keywords internal
.plot_bridge_outliers <- function(x) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  is_out  <- x$is_outlier & !is.na(x$is_outlier)
  ep_para <- x$bridge_eta_para
  ep_perp <- x$bridge_eta_perp

  has_dir <- !is.null(ep_para) && !is.null(ep_perp) &&
             any(is.finite(ep_para)) && any(is.finite(ep_perp))

  ## left: directional scatter when available; sorted-eta fallback otherwise
  if (has_dir) {
    finite <- is.finite(ep_para) & is.finite(ep_perp)
    xs <- log10(pmax(ep_para[finite], 0) + 1)
    ys <- log10(pmax(ep_perp[finite], 0) + 1)
    flag_idx <- finite & is_out
    graphics::plot(xs, ys,
                   xlab = expression(log[10](eta["para"] + 1)),
                   ylab = expression(log[10](eta["perp"] + 1)),
                   pch = 16, cex = 0.4,
                   col = grDevices::adjustcolor("grey50", 0.5),
                   main = "Bridge directional decomposition")
    graphics::abline(0, 1, lty = 2, col = "grey40")
    if (any(flag_idx)) {
      graphics::points(log10(pmax(ep_para[flag_idx], 0) + 1),
                       log10(pmax(ep_perp[flag_idx], 0) + 1),
                       pch = 1, cex = 1.4, col = "firebrick", lwd = 1.4)
    }
  } else {
    eta <- x$bridge_eta
    le <- log(eta[is.finite(eta) & eta > 0])
    if (length(le) >= 10) {
      ol <- sort(le)
      graphics::plot(ol, seq_along(ol), type = "l",
                     xlab = expression(log(eta)),
                     ylab = "sorted rank",
                     main = "Sorted bridge scores")
      if (any(is_out, na.rm = TRUE)) {
        out_le <- log(eta[is_out & is.finite(eta) & eta > 0])
        graphics::abline(v = max(out_le, na.rm = TRUE) + 0.01,
                         lty = 2, col = "red")
      }
    } else {
      graphics::plot.new()
    }
  }

  ## right: map with flagged points
  cc <- sf::st_coordinates(x)
  graphics::plot(cc, type = "l", col = "grey70", asp = 1,
                 xlab = "x (m)", ylab = "y (m)",
                 main = sprintf("Bridge flags (%d)", sum(is_out)))
  graphics::points(cc, col = "grey40", cex = 0.2, pch = 19)
  if (any(is_out, na.rm = TRUE)) {
    graphics::points(cc[is_out, , drop = FALSE],
                     col = "red", pch = 4, cex = 0.9, lwd = 1.5)
  }
  invisible(NULL)
}

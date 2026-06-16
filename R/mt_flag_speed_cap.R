#' Flag fixes by a speed cap (data-driven or physiological)
#'
#' Flags fixes that straddle an implausibly fast step.  The threshold
#' can be derived data-driven from the track's own step-speed
#' distribution (entropy valley or dip-test-validated broken-stick
#' break) or supplied as a hard physiological cap in m/s.  This is
#' the only primitive that scores steps as first-class objects, and
#' therefore the only one that catches the boundary transitions of a
#' coherent outlier block that per-fix detectors structurally cannot
#' see.
#'
#' @details
#' For each fix \eqn{i}, let \eqn{v_i = d_i / \Delta t_i} be the
#' implied step speed, where \eqn{d_i} is the distance to fix
#' \eqn{i+1} (via \code{\link[move2]{mt_distance}}, geodesic for
#' lon/lat data, Euclidean for projected data) and \eqn{\Delta t_i}
#' the corresponding time lag.  A fix is flagged
#' \code{is_speed_above_cap = TRUE} if either its outgoing step
#' \eqn{v_i} or its incoming step \eqn{v_{i-1}} exceeds the cap ---
#' without more information the detector cannot tell which of the
#' two fixes the offending geometry belongs to, so both endpoints
#' are flagged.
#'
#' The cap itself is chosen by \code{threshold_type}:
#' \itemize{
#'   \item \code{"auto"} (default): entropy-valley on \eqn{-\log v};
#'         if none is found, fall back to the broken-stick / gap
#'         break, retained only if Hartigan's dip test confirms the
#'         distribution is significantly multimodal (\eqn{p < 0.05}).
#'         On clean unimodal speed distributions this returns no
#'         flags --- the honest "no outliers" outcome.  The most
#'         principled choice, and the default.
#'   \item \code{"entropy"}: entropy valley only.  Most conservative.
#'   \item \code{"gap"}: broken-stick plus tail-decay inflection only.
#'         More sensitive but can over-flag on legitimate heavy tails.
#'   \item \code{"hard"}: uses \code{v_max} as a literal physiological
#'         cap.  Requires \code{v_max} to be supplied.  Pick species'
#'         documented biomechanical top speed plus ~25\% margin.
#' }
#'
#' When \code{jitter} is supplied, an analogous \emph{lower} bound
#' is applied to step length (not speed): \code{is_step_below_jitter
#' = TRUE} where either neighbouring step is below \code{jitter}
#' metres.  Use this only if you explicitly want to mark sub-noise
#' floor displacements.
#'
#' On multi-track input the function dispatches per individual: each
#' track's data-driven cap is tuned to that track's own speed
#' distribution, and \code{attr(out, "v_max_used")} is returned as a
#' named numeric (one entry per track id).  Single-track input retains
#' the scalar attribute.  Empty geometries and non-finite (\code{NA})
#' coordinates are a hard error (remove them upstream); rows with missing
#' timestamps yield NA step speeds and are excluded from flagging but
#' retained in the output, preserving row parity with the caller.
#'
#' @param x A \code{move2} object.  Longitude/latitude or projected;
#'   step lengths are computed via the package's geodesic / Euclidean
#'   helper which handles both CRS types directly.
#' @param v_max Numeric scalar in m/s or \code{NULL} (default).
#'   Required when \code{threshold_type = "hard"}; ignored otherwise.
#'   When \code{threshold_type != "hard"} the cap is computed from
#'   the data.
#' @param threshold_type Character, one of \code{"auto"} (default),
#'   \code{"entropy"}, \code{"gap"}, or \code{"hard"}.  See details.
#' @param threshold Numeric tuning parameter passed to the underlying
#'   threshold helper (entropy valley-depth ratio or gap multiplier).
#'   \code{NULL} (default) uses that helper's own default.  Ignored
#'   when \code{threshold_type = "hard"}.
#' @param jitter Numeric scalar or \code{NULL} (default).  Optional
#'   lower bound on absolute step length (metres).
#' @param physiological_ceiling Numeric scalar in m/s, or \code{NULL}
#'   (default).  Soft upper-bound used by the auto-cap path's
#'   biological-sanity warning: when the data-driven cap exceeds
#'   this value, a message is emitted suggesting the user supply a
#'   hard \code{v_max} or the \code{(mass, mode)} allometric prior
#'   to \code{mt_clean_track()}.  \code{NULL} falls back to the
#'   universal mode-agnostic ceiling of 55 m/s (Hirt et al. 2017
#'   95\% upper-CI of the maximum biological speed across all masses
#'   and modes).  When the user has \code{(mass, mode)} information,
#'   passing \code{physiological_ceiling = v_phys_estimate(mass,
#'   mode) * 1.25} gives a sharper per-species check (sprint margin
#'   on the central allometric prediction).  The check is
#'   warning-only; it never alters the cap.  Ignored when
#'   \code{threshold_type = "hard"} (the user's cap stands).
#' @param pool_by Optional character vector of length 1 or 2 naming
#'   column(s) in \code{mt_track_data(x)}.  Length 1: single column
#'   used as both fit set and operating unit.  Length 2:
#'   \code{c(outer, inner)} where \code{outer} names the fit-source
#'   column (the union of its events supplies the step-speed
#'   distribution that \code{.compute_speed_cap} -- with its
#'   fraction-above and mode-position gates -- runs over) and
#'   \code{inner} names the operating unit (within which the
#'   pool-fitted cap is applied and flags are unioned).  Length 2
#'   requires strict nesting: every distinct \code{inner} value
#'   must map to exactly one \code{outer} value.  Length \eqn{> 2}
#'   is rejected.  Pool flags union into \code{is_outlier} +
#'   \code{is_speed_above_cap} -- additive, never un-flags what
#'   per-track caught.  Per-track pool caps are exposed via
#'   \code{attr(out, "v_max_used_pool")}, alongside the per-track
#'   \code{attr(out, "v_max_used")}.  \code{NULL} (default)
#'   preserves per-track behaviour byte-identically.  Ignored
#'   when \code{threshold_type = "hard"} (the user-supplied cap
#'   is already uniform across tracks).
#' @param plot Logical.  If \code{TRUE} (default), produce a diagnostic
#'   map with flagged fixes highlighted.
#' @param remove Logical.  If \code{TRUE}, drop flagged rows from the
#'   returned object.  Default \code{FALSE}.
#' @param silent Logical.  If \code{FALSE} (default) the function
#'   prints a one-line summary of the cap chosen, gate decisions, and
#'   the count of flagged fixes.  Set \code{TRUE} to suppress.  Errors
#'   and warnings are always shown.
#'
#' @return The input object with added columns:
#'   \describe{
#'     \item{\code{step_speed}}{Implied outgoing step speed in m/s.
#'           \code{NA} for the last fix of each track.}
#'     \item{\code{is_speed_above_cap}}{Logical.  \code{TRUE} where
#'           either adjacent step exceeds the cap.}
#'     \item{\code{is_step_below_jitter}}{Logical.  Present only when
#'           \code{jitter} is supplied.}
#'     \item{\code{is_outlier}}{Logical.  Union of the above flags.}
#'     \item{\code{flagged_by_speed}}{Logical.  Same as
#'       \code{is_outlier} for this primitive; named for parity with the
#'       other primitives' output schema (so the output can be voted by
#'       \code{\link{mt_flag_consensus}}).}
#'     \item{\code{loglr_speed}}{Signed log-likelihood-ratio of outlier
#'           vs not for the above-cap signal, in nat units
#'           (\code{log(step_speed / cap)}); \code{> 0} above the cap.
#'           Two-sided in the data-driven modes; one-sided (negative
#'           clamped to 0) under \code{threshold_type = "hard"}, where a
#'           slow step is no evidence.  See
#'           \code{DESIGN_evidence_accumulation.md}.}
#'   }
#'   The value of the cap actually used is stored as the attribute
#'   \code{"v_max_used"}.
#'
#' @examples
#' \dontrun{
#' library(move2)
#' ## data-driven (default); silent on unimodal speed distributions
#' x <- mt_flag_speed_cap(x)
#'
#' ## hard physiological cap: golden eagle documented top ~50 m/s
#' x <- mt_flag_speed_cap(x, v_max = 50, threshold_type = "hard")
#'
#' ## diagnostic: inspect distribution + suggested cap without flagging
#' v <- mt_suggest_speed_cap(x)
#' }
#'
#' @seealso \code{\link{mt_suggest_speed_cap}} for a diagnostic-only
#'   helper that returns the suggested value without modifying
#'   \code{x}; \code{\link{mt_flag_outliers_bridge}} and
#'   \code{\link{mt_flag_outliers}} for per-fix detectors.
#'
#' @importFrom move2 mt_distance mt_time_lags mt_track_id
#' @importFrom graphics par plot points legend
#' @importFrom grDevices adjustcolor
#' @export
mt_flag_speed_cap <- function(x,
                               v_max                  = NULL,
                               threshold_type         = c("auto", "entropy",
                                                          "gap", "hard"),
                               threshold              = NULL,
                               jitter                 = NULL,
                               physiological_ceiling  = NULL,
                               pool_by                = NULL,
                               plot                   = TRUE,
                               remove                 = FALSE,
                               silent                 = FALSE) {

  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  threshold_type <- match.arg(threshold_type)
  if (!is.null(pool_by)) {
    ## `.resolve_pool_groups` validates shape, columns, and (for
    ## length-2 input) the strict-nesting requirement.  Validate
    ## unconditionally; bad inputs must error even when
    ## threshold_type = "hard" short-circuits the pool path.
    invisible(.resolve_pool_groups(x, pool_by, silent = TRUE))
  }
  ## brief narrator helper -- suppressed under silent = TRUE
  say <- function(...) if (!silent) message(...)
  ## physiological ceiling validation: NULL (use 55 m/s fallback) or
  ## a positive scalar.  The 55 m/s fallback is the Hirt 2017 mode-
  ## agnostic upper-CI; users with (mass, mode) should pass a
  ## sharper value (v_phys_estimate(mass, mode) * 1.25 is the
  ## suggested form).
  if (!is.null(physiological_ceiling)) {
    if (!is.numeric(physiological_ceiling) ||
        length(physiological_ceiling) != 1L ||
        is.na(physiological_ceiling) || physiological_ceiling <= 0) {
      rlang::abort(
        "`physiological_ceiling` must be NULL or a positive scalar (m/s).",
        class = "move2utils_mt_flag_speed_cap_bad_physiological_ceiling")
    }
  }
  if (threshold_type == "hard") {
    if (is.null(v_max) || !is.numeric(v_max) || length(v_max) != 1 ||
        is.na(v_max) || v_max <= 0) {
      rlang::abort(
        "`v_max` must be a positive scalar when threshold_type = \"hard\".",
        class = "move2utils_mt_flag_speed_cap_bad_v_max_for_hard")
    }
  } else if (!is.null(v_max)) {
    rlang::abort(paste0(
      "`v_max` is only used when threshold_type = \"hard\". ",
      "Either set threshold_type = \"hard\" or leave v_max = NULL."),
      class = "move2utils_mt_flag_speed_cap_v_max_without_hard")
  }
  if (!is.null(jitter)) {
    if (!is.numeric(jitter) || length(jitter) != 1 ||
        is.na(jitter) || jitter < 0) {
      rlang::abort("`jitter` must be a non-negative scalar (metres) or NULL.",
                   class = "move2utils_mt_flag_speed_cap_bad_jitter")
    }
  }

  ## ---- multi-track dispatch -------------------------------------
  ## For data-driven thresholds the cap should be tuned to each track's
  ## own speed distribution; pooling obscures per-individual scale.
  ## Per-track caps are returned as a named numeric on the merged
  ## object's "v_max_used" attribute.  Multi-track handling stays at
  ## wrapper level because the cap-merge logic is custom (not a
  ## generic post-rbind operation the dispatcher handles).
  ids <- move2::mt_track_id(x)
  unique_ids <- unique(ids)
  if (length(unique_ids) > 1L) {
    say("Processing ", length(unique_ids), " individuals separately...")
    results <- lapply(unique_ids, function(id) {
      xi <- x[ids == id, ]
      say(sprintf("--- %s (%d locations) ---", id, nrow(xi)))
      mt_flag_speed_cap(xi,
                        v_max                 = v_max,
                        threshold_type        = threshold_type,
                        threshold             = threshold,
                        jitter                = jitter,
                        physiological_ceiling = physiological_ceiling,
                        plot                  = FALSE,
                        remove                = FALSE,
                        silent                = silent)
    })
    out <- do.call(rbind, results)
    per_track_caps <- vapply(results,
                              function(r) {
                                v <- attr(r, "v_max_used")
                                if (is.null(v)) NA_real_ else as.numeric(v)
                              },
                              numeric(1))
    names(per_track_caps) <- as.character(unique_ids)
    attr(out, "v_max_used")     <- per_track_caps
    attr(out, "threshold_type") <- threshold_type

    ## ---- pool_by union: refit cap on union of step speeds per group ----
    ## Only runs in data-driven mode (auto / entropy / gap).  In "hard"
    ## mode the user-supplied v_max is uniform across tracks; pooling
    ## does nothing.  Pool flags union into is_outlier +
    ## is_speed_above_cap -- additive, never un-flag what per-track caught.
    if (!is.null(pool_by) && threshold_type != "hard") {
      out <- .speed_cap_pool_union(out, x, pool_by,
                                    threshold_type = threshold_type,
                                    threshold      = threshold,
                                    silent         = silent)
    }

    if (plot) .plot_speed_cap(out, stats::median(per_track_caps, na.rm = TRUE))
    out <- .attach_speedcap_loglr(out, threshold_type)
    if (remove) out <- out[!out$is_outlier, ]
    return(out)
  }

  was_longlat <- isTRUE(sf::st_is_longlat(x))

  ## ---- common plumbing: single-track + extraction + .speed_cap_fn_core
  ## + lift.  Stays in input CRS (Haversine on lon/lat, Euclidean on
  ## projected) for large-track performance. ----
  out <- .clean_track_dispatch(
    x,
    fn_core             = .speed_cap_fn_core,
    fn_core_args        = list(was_longlat          = was_longlat,
                                v_max                 = v_max,
                                threshold_type        = threshold_type,
                                threshold             = threshold,
                                jitter                = jitter,
                                physiological_ceiling = physiological_ceiling,
                                silent                = silent),
    per_track_args      = list(),
    need_time           = TRUE,
    hygiene_strict_time = FALSE,
    project_longlat     = FALSE,
    n_min               = 2L,
    n_min_severity      = "say",
    primitive_label     = "speed_cap",
    silent              = silent)

  if (plot) {
    v_used <- attr(out, "v_max_used")
    .plot_speed_cap(out, if (is.null(v_used)) NA else v_used)
  }
  out <- .attach_speedcap_loglr(out, threshold_type)
  if (remove) out <- out[!out$is_outlier, ]
  out
}

## Speed-cap log-LR, based on the above-cap signal (the primary
## physiological / density speed outlier; the jitter flag is a separate
## low-speed signal, not represented here).  Surprisal = log(step_speed)
## nats, boundary at the cap (least-fast flagged step) -> loglr =
## log(step_speed / cap).  In the data-driven modes this is two-sided (a
## step well below the cap is evidence it is fine).  In `hard` mode it is
## ONE-sided: a step below a physiological cap says nothing (many outliers
## are slow), so negative evidence is clamped to 0.
## @keywords internal
.attach_speedcap_loglr <- function(obj, threshold_type) {
  ## flag-source column for votability (== is_outlier for a single detector).
  obj$flagged_by_speed <- obj$is_outlier
  s    <- obj$step_speed
  flag <- obj$is_speed_above_cap
  if (is.null(s) || is.null(flag)) return(obj)
  neglogp <- log(s)
  neglogp[!is.finite(neglogp) | !is.finite(s) | s <= 0] <- NA_real_
  lr <- .loglr_grouped(neglogp, flag, move2::mt_track_id(obj))
  if (identical(threshold_type, "hard")) lr <- pmax(lr, 0)
  obj$loglr_speed <- lr
  obj
}


## Raw-matrix entry point used by mt_clean_track to skip per-iter
## sf-class slicing.  Operates on (cc, t_s, active_idx, was_longlat,
## ...) and returns active-indexed is_outlier + step_speed +
## is_speed_above_cap (+ is_step_below_jitter if jitter > 0) plus
## scalar attributes v_max_used + threshold_type.
##
## @keywords internal
.speed_cap_fn_core <- function(cc, t_s, active_idx, was_longlat,
                                v_max                  = NULL,
                                threshold_type         = "auto",
                                threshold              = NULL,
                                jitter                 = NULL,
                                physiological_ceiling  = NULL,
                                silent                 = TRUE) {

  ## Contract: active_idx sorted ascending with no duplicates.
  stopifnot(!is.unsorted(active_idx), !anyDuplicated(active_idx))

  say <- function(...) if (!silent) message(...)
  n_a  <- length(active_idx)
  cc_a <- cc[active_idx, , drop = FALSE]
  t_a  <- t_s[active_idx]

  ## step lengths + time lags on the active subset
  step_m <- .step_lengths_from_cc(cc_a, was_longlat)
  dt_s   <- c(diff(t_a), NA_real_)

  speed_out <- ifelse(is.finite(step_m) & is.finite(dt_s) & dt_s > 0,
                      step_m / dt_s, NA_real_)
  speed_in  <- c(NA_real_, utils::head(speed_out, -1L))
  step_in   <- c(NA_real_, utils::head(step_m,    -1L))

  ## ---- compute the cap ----
  v_max_used <- if (threshold_type == "hard") {
    v_max
  } else {
    .compute_speed_cap(speed_out, threshold_type, threshold, silent = silent)
  }

  ## ---- biological-sanity ceiling on data-driven caps ----
  ## Same defence-in-depth mirroring as before -- inner warning gets
  ## suppressed when mt_clean_track calls this via suppressMessages,
  ## the orchestrator emits the outer warning at iter == 1L.
  ceiling_used   <- if (is.null(physiological_ceiling))
                      55 else physiological_ceiling
  ceiling_source <- if (is.null(physiological_ceiling))
                      "Hirt 2017 universal upper-CI (~52.6 m/s, fastest flier)"
                    else
                      "user-supplied physiological_ceiling"
  if (threshold_type != "hard" && is.finite(v_max_used) &&
      v_max_used > ceiling_used) {
    say(sprintf(paste0(
      "Auto-cap landed at %.1f m/s -- above %.1f m/s (%s).  The gap ",
      "finder is detecting a structural break within the outlier tail ",
      "rather than between bulk and outliers.  Supply `(mass, mode)` ",
      "to `mt_clean_track()` (or pass a hard `v_max`) for a principled ",
      "physiological cap.  See `?v_phys_estimate`."),
      v_max_used, ceiling_used, ceiling_source))
  }

  ## ---- flags ----
  above_cap <- if (is.finite(v_max_used)) {
    (is.finite(speed_out) & speed_out > v_max_used) |
      (is.finite(speed_in)  & speed_in  > v_max_used)
  } else {
    rep(FALSE, n_a)
  }

  has_jitter <- !is.null(jitter)
  below_jitter <- if (has_jitter) {
    (is.finite(step_m)  & step_m  < jitter) |
      (is.finite(step_in) & step_in < jitter)
  } else {
    rep(FALSE, n_a)
  }

  is_outlier <- above_cap | below_jitter

  ## ---- report ----
  n_cap <- sum(above_cap, na.rm = TRUE)
  cap_desc <- if (is.finite(v_max_used)) {
    sprintf("%g m/s (%s)", v_max_used, threshold_type)
  } else {
    sprintf("no cap (threshold_type=%s, no structural break found)",
            threshold_type)
  }
  msg <- sprintf("Speed cap: %s -- %d fix(es) flagged", cap_desc, n_cap)
  if (has_jitter) {
    n_jit <- sum(below_jitter, na.rm = TRUE)
    msg <- sprintf("%s (+ %d jitter below %g m)", msg, n_jit, jitter)
  }
  msg <- sprintf("%s.  Total is_outlier = %d (%.3f%%).",
                 msg, sum(is_outlier), 100 * sum(is_outlier) / n_a)
  say(msg)

  out_list <- list(
    is_outlier         = is_outlier,
    step_speed         = speed_out,
    is_speed_above_cap = above_cap,
    v_max_used         = v_max_used,
    threshold_type     = threshold_type
  )
  if (has_jitter) {
    out_list$is_step_below_jitter <- below_jitter
  }
  out_list
}


## Compute a data-driven speed cap via the structural-break machinery.
## Applied to -log(speed) so low values correspond to high (outlier)
## speeds, matching the .entropy_threshold_lower / .gap_threshold_lower
## convention used across the package.
##
## Returns the suggested v_max in m/s, or Inf if no structural break
## is detected (i.e., the distribution is unimodal -- safe on clean
## data, nothing gets flagged).
##
## Two gates protect against the cap landing inside a legitimate
## activity mode rather than above the bulk:
##   1. fraction-above: the cap must isolate <= max_flag_fraction of
##      fixes (existing safety, prevents runaway on heavily-contaminated
##      data).
##   2. mode-position: the cap must lie strictly above the rightmost
##      substantive mode of log(speed) -- otherwise it is severing
##      a real activity mode rather than an outlier tail.  See
##      .gate_speed_cap_mode_position for the full rationale.
##
## Mode-position gate for the auto speed cap.
##
## Asks the question the fraction-above guard cannot: is the proposed
## cap separating an outlier tail from the bulk, or is it cutting
## into a legitimate activity mode (e.g. flight)?
##
## The principled disambiguation is the position of the cap relative
## to the rightmost SUBSTANTIVE mode of the log-speed distribution:
##
##   - A "substantive" mode is one whose basin of attraction (region
##     between the two flanking density valleys) contains at least
##     `mode_min_frac` of the data (default 2%).  Sparse modes
##     created by a handful of outliers (e.g. spoofs at 200 m/s on a
##     stork track, ~0.05% of fixes per spoof mode) are not
##     substantive and are correctly ignored, while small but real
##     activity modes (e.g. a stork that flies for ~3% of the day)
##     are correctly recognised as part of the bulk.  The 2% choice
##     is a tradeoff: tracks where 2--5% of fixes are *genuine*
##     outliers are unlikely candidates for auto-cap detection
##     anyway and should use a hard physiological cap.
##
##   - The rightmost substantive mode is the upper end of the bulk
##     distribution -- the highest speed regime the animal actually
##     reaches in normal behaviour.
##
##   - A cap above that mode is isolating an outlier tail.  A cap at
##     or below that mode is severing legitimate movement.
##
## On the synthetic CPF_A track (where the gate must allow), the
## rightmost substantive mode is around 21 m/s (real flight) and the
## auto cap lands at ~60 m/s -- the gate passes.  On the white-stork
## Mia track (where the gate must decline), the rightmost substantive
## mode is at ~3.5 m/s (the bird's flight regime) and the auto cap
## drifts down to 2.3 m/s -- the gate refuses.
##
## @keywords internal
.gate_speed_cap_mode_position <- function(speed, v_b,
                                            mode_min_frac    = 0.02,
                                            mode_min_density = 0.05,
                                            n_grid = 512,
                                            margin_bw = 0.5) {
  s <- speed[is.finite(speed) & speed > 0]
  n <- length(s)
  if (n < 30L) {
    return(list(allow = TRUE,
                reason = "too few positive speeds for mode detection",
                rightmost_mode = NA_real_, position = NA_character_))
  }
  log_s <- log(s)
  d <- stats::density(log_s, n = n_grid)
  bw <- d$bw

  ## local maxima and minima on the density curve
  dy        <- diff(d$y)
  is_max    <- which(c(FALSE, dy[-length(dy)] > 0 & dy[-1L] < 0, FALSE))
  is_min    <- which(c(FALSE, dy[-length(dy)] < 0 & dy[-1L] > 0, FALSE))

  if (length(is_max) == 0L) {
    ## Pathological density (no interior maximum) -- can't reason
    ## about modes; fall back to allow.
    return(list(allow = TRUE,
                reason = "no interior modes detected",
                rightmost_mode = NA_real_, position = NA_character_))
  }
  global_peak <- max(d$y[is_max])

  ## For each candidate mode, the basin of attraction is bounded by
  ## the surrounding density minima (or the data extremes).  A mode
  ## is "substantive" only when BOTH (a) its basin contains at least
  ## mode_min_frac of the data and (b) its peak density is at least
  ## mode_min_density times the global peak density.  The basin-
  ## fraction criterion alone treats a tight cluster of repeated
  ## outliers (e.g. three identical spike speeds at 83 m/s) as a
  ## "mode"; the density-ratio criterion correctly identifies them
  ## as low-density spikes vs the high-density bulk and rejects.
  basin_bounds <- function(idx) {
    left  <- is_min[is_min < idx]
    right <- is_min[is_min > idx]
    lo <- if (length(left))  d$x[max(left)]   else -Inf
    hi <- if (length(right)) d$x[min(right)]  else  Inf
    c(lo, hi)
  }

  substantive_modes <- numeric(0)
  for (idx in is_max) {
    bb <- basin_bounds(idx)
    frac <- mean(log_s >= bb[1] & log_s <= bb[2])
    density_ratio <- d$y[idx] / global_peak
    if (frac >= mode_min_frac && density_ratio >= mode_min_density) {
      substantive_modes <- c(substantive_modes, exp(d$x[idx]))
    }
  }

  if (length(substantive_modes) == 0L) {
    ## No substantive mode detected -- distribution is too dispersed
    ## or pathological.  Conservatively allow (the fraction-above
    ## guard upstream is the safety net).
    return(list(allow = TRUE,
                reason = "no substantive mode found",
                rightmost_mode = NA_real_, position = NA_character_))
  }

  rightmost <- max(substantive_modes)
  ## Require the cap to sit at least margin_bw bandwidths above the
  ## rightmost substantive mode, on the log scale.
  threshold <- exp(log(rightmost) + margin_bw * bw)

  if (v_b > threshold) {
    return(list(allow = TRUE,
                reason = sprintf(
                  "cap %.2f m/s sits %.2f bandwidths above rightmost mode at %.2f m/s",
                  v_b, (log(v_b) - log(rightmost)) / bw, rightmost),
                rightmost_mode = rightmost, position = "above"))
  }

  ## Cap at or below the rightmost substantive mode -- it is severing
  ## a real activity mode, not isolating an outlier tail.
  position <- if (v_b < rightmost) "below" else "at"
  list(allow = FALSE,
       reason = sprintf(
         "cap %.2f m/s sits %s rightmost substantive mode at %.2f m/s",
         v_b, position, rightmost),
       rightmost_mode = rightmost, position = position)
}


## Pool-union closure for mt_flag_speed_cap.
##
## After multi-track dispatch produces a per-track step_speed column +
## per-track caps, this helper fits one cap per pool group from the
## union of step_speeds and reflags each track using the pooled cap.
## The fraction-above guard and mode-position gate inside
## .compute_speed_cap apply to the pooled distribution -- that is the
## correct semantics (the group's flag-fraction and the group's mode
## structure are what gate the pooled cap).
##
## Flag union: a fix flags pooled-above-cap iff EITHER its outgoing
## step speed (step_speed) OR its incoming step speed (= prior fix's
## step_speed within the same track) exceeds the pool cap.  Within-
## track lag is taken from the rbind'd output, which preserves
## within-track row order.
##
## Stores per-group pooled caps on attr(out, "v_max_used_pool").
##
## @keywords internal
.speed_cap_pool_union <- function(out, x, pool_by,
                                   threshold_type, threshold,
                                   silent = FALSE) {

  ids_event_all <- as.character(move2::mt_track_id(out))
  step_speed_all <- out$step_speed

  pool_step <- function(fit_idx, apply_idx) {
    ## Fit pool cap from OUTER group's union of valid step speeds.
    speeds_fit <- step_speed_all[fit_idx]
    cap_pool <- .compute_speed_cap(speeds_fit, threshold_type, threshold,
                                    silent = silent)
    if (!is.finite(cap_pool)) return(rep(FALSE, length(apply_idx)))

    ## Apply to INNER group.  For each track in the inner group,
    ## reconstruct speed_in (= lag of step_speed within track) and
    ## union flags.
    speeds_app <- step_speed_all[apply_idx]
    ids_g <- ids_event_all[apply_idx]
    above_pool <- logical(length(apply_idx))
    for (tid in unique(ids_g)) {
      rel <- which(ids_g == tid)
      sp_out <- speeds_app[rel]
      sp_in  <- c(NA_real_, sp_out[-length(rel)])
      above_pool[rel] <- (is.finite(sp_out) & sp_out > cap_pool) |
                          (is.finite(sp_in)  & sp_in  > cap_pool)
    }
    if (!silent) {
      n_new <- sum(above_pool & !out$is_outlier[apply_idx], na.rm = TRUE)
      message(sprintf(
        "  pool_by[%s]: cap = %.3g m/s, %d new fix(es) flagged across %d event(s).",
        paste(pool_by, collapse = ","), cap_pool, n_new, length(apply_idx)))
    }
    ## Stash the pool cap per track for the v_max_used_pool attribute
    ## after .apply_pool_union returns.  Each track in this inner
    ## group records the same cap (they all nest in the same outer).
    for (tid in unique(ids_g)) attr_pool[[tid]] <<- cap_pool
    above_pool
  }

  ## Build a per-track v_max_pool map (filled as pool_step runs).
  ## Multiple tracks in one group write the same value; harmless.
  attr_pool <- list()

  out <- .apply_pool_union(out, x, pool_by, pool_step,
                            flag_cols = c("is_outlier",
                                          "is_speed_above_cap"),
                            silent = silent)

  ## Attach per-track pool cap as an additional attribute (does not
  ## replace v_max_used, which remains the per-track cap).
  if (length(attr_pool)) {
    pc_named <- unlist(attr_pool, use.names = TRUE)
    attr(out, "v_max_used_pool") <- pc_named
  }
  out
}


## @keywords internal
##
## max_flag_fraction default 0.05 is HEURISTIC -- the safety guard
## rejects a candidate cap that would flag >5% of the distribution on
## the principle that "if a single threshold would flag >5% of
## fixes, it is almost certainly cutting between bulk activity modes
## (rest vs movement) rather than separating outliers."  Plausible
## range: 0.02--0.10.  Lower values are more conservative.  Same
## structural reason as max_flag_fraction in mt_clean_track but
## different default because this gate operates per-call without
## iteration recovery.
.compute_speed_cap <- function(speed, threshold_type, threshold,
                                 max_flag_fraction = 0.05,
                                 silent = FALSE) {
  say <- function(...) if (!silent) message(...)
  s <- speed[is.finite(speed) & speed > 0]
  if (length(s) < 10L) return(Inf)
  nl <- -log(s)

  ## NULL `threshold` defers to the leaf formals (single source of
  ## truth).  Pre-2026-05-25 the auto branch forwarded literal 0.3 / 3
  ## explicitly to .entropy_or_dip_gap_threshold_lower, shadowing the
  ## leaf and breaking sweep override via assignInNamespace.  See
  ## audits/2026-05-25-parameter-propagation/findings.md §1.3.
  br <- switch(threshold_type,
    auto    = .entropy_or_dip_gap_threshold_lower(
      nl,
      entropy_threshold = threshold,  # NULL -> leaf default
      gap_threshold     = NULL,       # leaf default
      dip_alpha         = 0.05),
    entropy = if (is.null(threshold))
                .entropy_threshold_lower(nl)
              else
                .entropy_threshold_lower(nl, threshold = threshold),
    gap     = if (is.null(threshold))
                .gap_threshold_lower(nl)
              else
                .gap_threshold_lower(nl, threshold = threshold))

  if (is.na(br$break_value)) return(Inf)
  v <- exp(-br$break_value)

  ## ---- gate 1: fraction-above safety guard ----------------------
  frac_above <- mean(s > v, na.rm = TRUE)
  if (frac_above > max_flag_fraction) {
    say(sprintf(
      "Speed-cap auto gate: the data-driven break at %.2f m/s would flag %.1f%% of fixes -- this almost certainly separates bulk activity modes (rest vs movement), not outliers.",
      v, 100 * frac_above))
    say("  For heavily-contaminated data (>5% of fixes are real outliers), a hard physiological cap from species biology is the more reliable tool:")
    say(sprintf(
      "      mt_flag_speed_cap(x, v_max = <species top speed, m/s>, threshold_type = \"hard\")"))
    say("  Proceeding with no auto cap (v_max = Inf).")
    return(Inf)
  }

  ## ---- gate 2: mode-position guard ------------------------------
  ## A cap that sits at or below the rightmost substantive mode of
  ## the speed distribution is severing a real activity mode (e.g.
  ## flight) rather than isolating an outlier tail.  This catches
  ## the failure mode that the fraction-above guard misses: tracks
  ## where 3-4% of fixes are LEGITIMATELY in the upper mode (a real
  ## bird flying for a few hours of the day), so the cap fraction
  ## passes the 5% guard but lands inside that real mode.
  mp <- .gate_speed_cap_mode_position(s, v)
  if (!mp$allow) {
    say(sprintf(
      "Speed-cap auto gate: the data-driven break at %.2f m/s sits %s the rightmost substantive mode at %.2f m/s -- the cap is separating activity modes (rest vs movement), not outliers from valid movement.",
      v, mp$position, mp$rightmost_mode))
    say("  For a fast-moving species, supply a physiological cap via:")
    say(sprintf(
      "      mt_flag_speed_cap(x, v_max = <species top speed, m/s>, threshold_type = \"hard\")"))
    say("  Or pass (mass, mode) to mt_clean_track for an allometric prior.  Proceeding with no auto cap (v_max = Inf).")
    return(Inf)
  }
  v
}


#' Suggest a data-driven speed cap from the observed step-speed distribution
#'
#' Companion diagnostic for \code{\link{mt_flag_speed_cap}}: reads the
#' distribution of implied step speeds on a track, looks for a structural
#' break separating plausible-speed fixes from extreme-speed outliers,
#' and suggests a value for \code{v_max}.  The suggestion is data-driven
#' (entropy-valley or broken-stick break in log-speed space) and
#' therefore still distribution-free; users can override with domain
#' knowledge if the suggested value disagrees with species biology.
#'
#' @details
#' K02-style spoof-boundary transitions produce speeds that sit in a
#' yawning gap above the bulk of the distribution: on that track,
#' 99.99\% of steps are below 32 m/s and the tail jumps directly to
#' 209 m/s with nothing in between.  Any reasonable break-detector
#' places the suggested cap in that gap.  When no such gap exists
#' (clean data without extreme-speed outliers), the function returns
#' \code{NA} and reports that no structural break was found --- the
#' honest "no outliers" outcome.
#'
#' \strong{Suggester / flagger alignment (v0.3).}  The suggester now
#' applies the same pipeline as \code{\link{mt_flag_speed_cap}}'s auto
#' path: identical default \code{method = "auto"}, identical break
#' detection, identical \code{>5\%} fraction-above guard, identical
#' mode-position gate (cap must sit above the rightmost substantive
#' mode of the speed distribution), identical 55 m/s biological-sanity
#' ceiling warning.  A user inspecting their data with
#' \code{mt_suggest_speed_cap()} and then calling
#' \code{mt_flag_speed_cap(x, threshold_type = "auto")} now sees a
#' single coherent recommendation rather than two independent
#' proposals.  Pre-v0.3 versions of this function used a different
#' default method and skipped the mode-position gate, leading to
#' silent disagreement between the suggester and the flagger on
#' multi-state tracks (the "rest vs flight" failure mode).
#'
#' @param x A \code{move2} object.
#' @param method Character, one of \code{"auto"} (default, matches
#'   \code{\link{mt_flag_speed_cap}}'s auto path: entropy-valley first,
#'   broken-stick + dip-test fallback), \code{"entropy"} (entropy-
#'   valley only; conservative), or \code{"gap"} (broken-stick + tail-
#'   decay only; sensitive).  Prior to v0.3 the default was
#'   \code{"entropy"}; the change to \code{"auto"} brings the suggester
#'   into agreement with the flagger so a user inspecting their data
#'   then calling \code{\link{mt_flag_speed_cap}} sees the same
#'   pipeline.
#' @param threshold Numeric passed to the underlying threshold helper.
#'   Default \code{NULL} uses that helper's own default.
#' @param physiological_ceiling Optional numeric in m/s.  When supplied
#'   (typically as \code{v_phys_estimate(mass, mode) * 1.25} -- a
#'   sprint-margined per-species ceiling), the function warns when the
#'   suggested cap exceeds this value, matching
#'   \code{\link{mt_flag_speed_cap}}'s warning text.  When \code{NULL},
#'   falls back to the universal mode-agnostic 55 m/s ceiling derived
#'   from Hirt 2017 (the upper-95\%-CI of the fastest biological
#'   flier).  A suggested cap above the ceiling indicates the gap
#'   finder is detecting a structural break inside the data's outlier
#'   tail rather than between bulk and outliers.
#' @param mass,mode Optional.  When \code{mass} is supplied (kg), an
#'   allometric prediction of the species' physiological maximum
#'   speed is derived via \code{\link{v_phys_estimate}} (Hirt et
#'   al. 2017) and overlaid on the diagnostic plot.  If \code{mode}
#'   is also given (one of \code{"flying"}, \code{"running"},
#'   \code{"swimming"}), a single mode-specific line is shown.  If
#'   \code{mode} is left \code{NULL}, all three mode-specific lines
#'   are shown so the user can compare predictions for ambiguous
#'   cases (mode-switching species, uncertain locomotor categorisation).
#' @param v_max Optional numeric scalar in m/s.  If supplied, marks
#'   the user's own physiological cap on the diagnostic plot --- the
#'   same value they would pass to \code{\link{mt_flag_speed_cap}} or
#'   \code{\link{mt_clean_track}}.  When both \code{(mass, mode)} and
#'   \code{v_max} are given, the diagnostic shows three independent
#'   perspectives on the same question (empirical, allometric, user)
#'   so the user can read their congruence directly.
#' @param plot Logical.  If \code{TRUE} (default), render a density
#'   diagnostic with the suggested cap overlaid as a vertical line,
#'   plus the quantile table printed to \code{message()}.
#' @param silent Logical.  If \code{FALSE} (default) the function
#'   prints quantile, congruence, and triangulation diagnostics.
#'   Set \code{TRUE} to suppress.
#'
#' @return A numeric scalar: the suggested \code{v_max} in m/s, or
#'   \code{NA_real_} if no structural break was found.  Returned
#'   invisibly when \code{plot = TRUE}.
#'
#' @examples
#' \dontrun{
#' v <- mt_suggest_speed_cap(track)
#' if (!is.na(v)) track <- mt_flag_speed_cap(track, v_max = v)
#'
#' ## triangulate empirical, allometric and user estimates:
#' mt_suggest_speed_cap(track, mass = 5, mode = "flying", v_max = 30)
#' }
#'
#' @importFrom stats density quantile setNames
#' @importFrom graphics abline axis hist par plot
#' @export
mt_suggest_speed_cap <- function(x, method = c("auto", "entropy", "gap"),
                                  threshold = NULL,
                                  mass = NULL, mode = NULL,
                                  v_max = NULL,
                                  physiological_ceiling = NULL,
                                  plot = TRUE,
                                  silent = FALSE) {
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  method <- match.arg(method)
  if (!is.null(physiological_ceiling)) {
    if (!is.numeric(physiological_ceiling) ||
        length(physiological_ceiling) != 1L ||
        is.na(physiological_ceiling) || physiological_ceiling <= 0) {
      rlang::abort(
        "`physiological_ceiling` must be a positive scalar (m/s) or NULL.",
        class = "move2utils_mt_suggest_speed_cap_bad_physiological_ceiling")
    }
  }
  ## brief narrator helper -- suppressed under silent = TRUE
  say <- function(...) if (!silent) message(...)

  ## ---- optional allometric prior --------------------------------
  ## v_allo is a *named* numeric vector: when mode is supplied, length 1
  ## with the matching name; when mode is NULL but mass is given, length
  ## 3 with one entry per mode.
  v_allo <- numeric(0)
  if (!is.null(mass)) {
    if (!is.null(mode)) {
      allo <- v_phys_estimate(mass = mass, mode = mode)
      v_allo <- setNames(as.numeric(allo), mode)
    } else {
      v_allo <- vapply(c("flying", "running", "swimming"),
                       function(m) suppressWarnings(
                         as.numeric(v_phys_estimate(mass = mass, mode = m))),
                       numeric(1))
    }
  } else if (!is.null(mode)) {
    rlang::abort(
      "`mode` was supplied but `mass` is NULL; provide both, or just `mass` to see all three mode predictions.",
      class = "move2utils_mt_suggest_speed_cap_mode_without_mass")
  }
  if (!is.null(v_max)) {
    if (!is.numeric(v_max) || length(v_max) != 1L ||
        is.na(v_max) || v_max <= 0) {
      rlang::abort("`v_max` must be a positive scalar (m/s) or NULL.",
                   class = "move2utils_mt_suggest_speed_cap_bad_v_max")
    }
  }

  step_m <- .step_lengths_fast(x)
  dt_s   <- as.numeric(move2::mt_time_lags(x, units = "secs"))
  speed  <- ifelse(is.finite(step_m) & is.finite(dt_s) & dt_s > 0,
                   step_m / dt_s, NA_real_)
  speed  <- speed[is.finite(speed) & speed > 0]

  if (length(speed) < 10) {
    say("Too few positive step speeds for break detection.")
    return(invisible(NA_real_))
  }

  ## Quantile summary
  qs <- stats::quantile(speed,
                         c(0.5, 0.9, 0.99, 0.999, 0.9999, 1),
                         na.rm = TRUE)
  say("Step-speed quantiles (m/s):")
  for (i in seq_along(qs)) {
    say(sprintf("  %6s : %7.2f", names(qs)[i], qs[i]))
  }

  ## Apply the same break-detection + gating pipeline the flagger uses
  ## (.compute_speed_cap handles threshold dispatch, the >5%-fraction
  ## guard, and the mode-position gate in a single call).  Returns Inf
  ## when no break is found or when either gate refuses; we translate
  ## Inf to NA for the suggester's documented "no actionable suggestion"
  ## semantic.
  v_suggest <- .compute_speed_cap(speed,
                                    threshold_type = method,
                                    threshold      = threshold,
                                    silent         = silent)
  if (!is.finite(v_suggest)) v_suggest <- NA_real_

  ## Congruence diagnostic: run the OTHER threshold method on -log(speed)
  ## without the gating, so the user sees the raw structural information
  ## available from each detector.  A factor-of-3+ disagreement between
  ## the two raw estimates indicates multi-scale tail structure worth
  ## inspecting before setting a final v_max.
  nl <- -log(speed)
  if (method == "auto") {
    ## "auto" is entropy-then-gap; report both raw values for
    ## transparency.
    ## Defer to leaf formals -- the package-wide single source of
    ## truth for the entropy / gap defaults.
    br_entropy <- .entropy_threshold_lower(nl)
    br_gap     <- .gap_threshold_lower(nl)
    v_entropy <- if (is.na(br_entropy$break_value)) NA_real_
                 else exp(-br_entropy$break_value)
    v_gap     <- if (is.na(br_gap$break_value)) NA_real_
                 else exp(-br_gap$break_value)
    if (!is.na(v_entropy) && !is.na(v_gap)) {
      disagreement <- max(v_entropy, v_gap) / min(v_entropy, v_gap)
      say(sprintf(
        "Congruence: entropy raw = %.2f m/s vs gap raw = %.2f m/s -- disagreement %.1f x (gates not applied).",
        v_entropy, v_gap, disagreement))
      if (disagreement > 3) {
        say("  Large disagreement indicates multi-scale tail structure (e.g. coherent outlier cluster plus extreme singletons).")
        say("  Inspect the speed ECDF and consider both candidates when choosing v_max.")
      }
    } else if (!is.na(v_entropy)) {
      say(sprintf("Congruence: entropy raw = %.2f m/s; gap found no break.",
                  v_entropy))
    } else if (!is.na(v_gap)) {
      say(sprintf("Congruence: gap raw = %.2f m/s; entropy found no break.",
                  v_gap))
    }
  } else {
    ## entropy or gap selected explicitly -- compare to the other.
    other_fn <- switch(method,
      entropy = function(v) .gap_threshold_lower(v),     # leaf default
      gap     = function(v) .entropy_threshold_lower(v)) # leaf default
    br_other <- other_fn(nl)
    v_other <- if (is.na(br_other$break_value)) NA_real_
                else exp(-br_other$break_value)
    other_method <- if (method == "entropy") "gap" else "entropy"
    if (!is.na(v_suggest) && !is.na(v_other)) {
      disagreement <- max(v_suggest, v_other) / min(v_suggest, v_other)
      say(sprintf(
        "Congruence: %s = %.2f m/s vs %s raw = %.2f m/s -- disagreement %.1f x.",
        method, v_suggest, other_method, v_other, disagreement))
      if (disagreement > 3) {
        say("  Large disagreement indicates multi-scale tail structure (e.g. coherent outlier cluster plus extreme singletons).")
        say("  Inspect the speed ECDF and consider both candidates when choosing v_max.")
      }
    } else if (!is.na(v_other)) {
      say(sprintf("Congruence: %s found no break but %s raw suggests %.2f m/s; methods disagree.",
                       method, other_method, v_other))
    }
  }

  ## Biological-sanity ceiling warning: mirror the flagger's behaviour
  ## (mt_flag_speed_cap.R:243-282).  When the suggested cap exceeds the
  ## ceiling, the gap finder is detecting a structural break inside the
  ## data's outlier tail rather than between bulk and outliers -- the
  ## user should consult species biology rather than trust the suggestion.
  if (!is.na(v_suggest)) {
    ceiling_used   <- if (is.null(physiological_ceiling))
                        55 else physiological_ceiling
    ceiling_source <- if (is.null(physiological_ceiling))
                        "Hirt 2017 universal upper-CI (~52.6 m/s, fastest flier)"
                      else
                        "user-supplied physiological_ceiling"
    if (v_suggest > ceiling_used) {
      say(sprintf(paste0(
        "Suggested cap %.1f m/s exceeds %.1f m/s (%s).  The gap finder ",
        "is detecting a structural break within the outlier tail rather ",
        "than between bulk and outliers.  Consult species biology (or ",
        "supply `(mass, mode)` to derive an allometric cap) rather than ",
        "trusting this suggestion."),
        v_suggest, ceiling_used, ceiling_source))
    }
  }

  if (is.na(v_suggest)) {
    say(sprintf("No structural break detected (method = %s). ",
                    method),
            "The distribution appears unimodal; no outlier-like tail. ",
            "If you still want a cap, set it from species biology.")
  } else {
    say(sprintf("Suggested v_max = %.2f m/s (method = %s, %d of %d ",
                    v_suggest, method,
                    sum(speed > v_suggest), length(speed)),
            sprintf("steps above cap, %.4f%%).",
                    100 * sum(speed > v_suggest) / length(speed)))
  }

  ## ---- triangulation diagnostic ----------------------------------
  ## When an allometric prior or user estimate is available, compare
  ## with the empirical break.  Disagreement is itself diagnostic.
  if (length(v_allo) || !is.null(v_max)) {
    say("Triangulation:")
    for (nm in names(v_allo)) {
      say(sprintf(
        "  allometric (Hirt 2017, %-8s mass=%g kg): %.2f m/s",
        nm, mass, v_allo[[nm]]))
    }
    if (!is.null(v_max)) {
      say(sprintf("  user-supplied v_max                        : %.2f m/s",
                      v_max))
    }
    if (!is.na(v_suggest)) {
      say(sprintf("  empirical break (%-8s)                  : %.2f m/s",
                      method, v_suggest))
    }
    ## When a single mode-specific allometric is given, run the
    ## congruence diagnostic against it.  When all three modes are
    ## shown, leave the comparison to the user since "the right mode"
    ## is itself their choice.
    if (length(v_allo) == 1L && !is.na(v_suggest)) {
      v_a <- as.numeric(v_allo)
      ratio <- v_suggest / v_a
      if (ratio > 2) {
        say(sprintf(
          "  Empirical break (%.2f m/s) lies %.1fx above the allometric prediction (%.2f m/s).",
          v_suggest, ratio, v_a))
        say("  This indicates contamination above the physiological cap; the allometric prediction is the more reliable cut.")
      } else if (ratio < 0.5) {
        say(sprintf(
          "  Empirical break (%.2f m/s) lies well below the allometric prediction (%.2f m/s).",
          v_suggest, v_a))
        say("  This indicates substantive within-distribution structure (e.g. behavioural state changes); cutting at the empirical break would over-flag.  The allometric prediction is the appropriate hard upper bound.")
      } else {
        say(sprintf(
          "  Empirical and allometric within %.1fx of each other -- cap is well-determined.",
          max(ratio, 1 / ratio)))
      }
    }
    if (!is.null(v_max) && length(v_allo) == 1L) {
      v_a <- as.numeric(v_allo)
      r2 <- v_max / v_a
      if (r2 > 2) {
        say(sprintf(
          "  User v_max (%.2f m/s) is %.1fx the allometric prediction.  Confirm intended (e.g. specialist sprint, peregrine stoop in 3D data).",
          v_max, r2))
      } else if (r2 < 0.5) {
        say(sprintf(
          "  User v_max (%.2f m/s) is more conservative than the allometric prediction (%.2f m/s).  Confirm intended.",
          v_max, v_a))
      }
    }
  }

  if (plot) {
    .plot_speed_suggest(speed, v_suggest, v_allo = v_allo, v_user = v_max)
  }

  invisible(v_suggest)
}


## Density plot of log-speed with vertical reference lines: the
## empirical break (red), one or three Hirt 2017 allometric
## predictions (blues -- one per locomotor mode), and the user's own
## v_max (green) when supplied.
##
## @keywords internal
.plot_speed_suggest <- function(speed, v_suggest,
                                  v_allo = numeric(0), v_user = NULL) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(4, 4, 3, 1))

  d <- stats::density(log(speed), bw = "nrd0", n = 512)
  graphics::plot(d, main = "Step-speed distribution",
                 xlab = "log(speed)  [log m/s]",
                 ylab = "density",
                 col = "grey30", lwd = 1.5)

  legends <- character()
  cols    <- character()
  ltys    <- integer()

  if (!is.na(v_suggest)) {
    graphics::abline(v = log(v_suggest), col = "red", lwd = 2, lty = 2)
    legends <- c(legends, sprintf("empirical break = %.1f m/s", v_suggest))
    cols    <- c(cols, "red")
    ltys    <- c(ltys, 2L)
  }
  ## allometric lines: one or three, distinct hues per mode
  allo_palette <- c(flying   = "steelblue4",
                    running  = "tan4",
                    swimming = "purple4")
  for (nm in names(v_allo)) {
    val <- v_allo[[nm]]
    if (is.finite(val) && val > 0) {
      graphics::abline(v = log(val), col = allo_palette[[nm]],
                       lwd = 2, lty = 4)
      legends <- c(legends, sprintf("Hirt %s = %.1f m/s", nm, val))
      cols    <- c(cols, allo_palette[[nm]])
      ltys    <- c(ltys, 4L)
    }
  }
  if (!is.null(v_user) && is.numeric(v_user) && !is.na(v_user) &&
      v_user > 0) {
    graphics::abline(v = log(v_user), col = "darkgreen",
                     lwd = 2, lty = 3)
    legends <- c(legends, sprintf("user v_max = %.1f m/s", v_user))
    cols    <- c(cols, "darkgreen")
    ltys    <- c(ltys, 3L)
  }
  if (length(legends)) {
    graphics::legend("topright", legend = legends, col = cols,
                     lty = ltys, lwd = 2, bty = "n", cex = 0.85)
  } else {
    graphics::legend("topright", legend = "no structural break found",
                     bty = "n")
  }
  invisible(NULL)
}


## Diagnostic plot: ECDF with the applied cap overlaid + track map
## with flagged fixes.
##
## Left panel: empirical CDF of step speeds on log-x, with the cap
## line and the n_above / frac_above annotation that makes the cost
## of the chosen cut visible.  Right panel: track outline with
## flagged fixes.
##
## @keywords internal
.plot_speed_cap <- function(x, v_max_used) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  speed <- x$step_speed
  speed <- speed[is.finite(speed) & speed > 0]
  if (length(speed) >= 10) {
    sp_sort <- sort(speed)
    ec <- seq_along(sp_sort) / length(sp_sort)
    graphics::plot(sp_sort, ec, type = "l", log = "x",
                   xlab = "step speed (m/s, log scale)",
                   ylab = "empirical CDF",
                   col = "grey30", lwd = 1.4,
                   main = "Step-speed ECDF")
    if (is.finite(v_max_used)) {
      n_above   <- sum(speed > v_max_used)
      frac_text <- sprintf("v_max = %g m/s\n%d above (%.2f%%)",
                            v_max_used, n_above,
                            100 * n_above / length(speed))
      graphics::abline(v = v_max_used, col = "firebrick",
                       lwd = 1.6, lty = 2)
      graphics::legend("bottomright", legend = frac_text,
                       col = "firebrick", lty = 2, lwd = 1.6, bty = "n")
    } else {
      graphics::legend("bottomright", legend = "no cap applied",
                       bty = "n")
    }
  } else {
    graphics::plot.new()
  }

  cc <- sf::st_coordinates(x)
  title_txt <- if (is.finite(v_max_used)) {
    sprintf("Speed-cap flags (v_max = %g m/s)", v_max_used)
  } else {
    "Speed-cap flags (no cap applied)"
  }
  graphics::plot(cc, type = "l", col = "grey70", lwd = 0.3,
                 xlab = "", ylab = "", main = title_txt, asp = 1)
  flagged <- which(x$is_outlier)
  if (length(flagged)) {
    graphics::points(cc[flagged, , drop = FALSE],
                     col = grDevices::adjustcolor("red", 0.6),
                     pch = 20, cex = 0.6)
    graphics::legend("topright",
                     legend = sprintf("flagged: %d", length(flagged)),
                     col = "red", pch = 20, bty = "n")
  }
  invisible(NULL)
}

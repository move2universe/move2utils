#' Combine per-detector flag columns into a single outlier decision
#'
#' Many outlier-detection workflows produce per-detector flags
#' (`flagged_by_bridge`, `flagged_by_prob`, `flagged_by_speed`,
#' `flagged_by_detour`) and then ask a separate question: given the
#' agreement structure across detectors, which fixes are outliers?
#' `mt_flag_consensus()` is the canonical home for that decision.
#'
#' The function is used internally by [mt_clean_track()] each iteration,
#' but is also exported so users running their own cascades (or running
#' the four primitives standalone) can apply the same consensus rules
#' to their own per-detector flag columns.
#'
#' @section Detector set:
#' The decision is taken over an extensible set of detectors named by
#' the `detectors` argument (a named character vector mapping a
#' detector name to the flag column that holds its decision). The
#' default names the four standard cascade detectors. Adding a detector
#' is adding an entry to `detectors`; the built-in modes that count
#' votes (`majority`, `any`, the `consensus` clause of `class_aware`)
#' generalise automatically, while the named class rules
#' (`geometric_spike`, etc.) reference detectors by name and ignore
#' columns they do not mention.
#'
#' @section Consensus modes:
#' \describe{
#'   \item{\code{"class_aware"}}{The default prior to v0.4.0
#'     (round-4 audit, 2026-05-11); the current default is
#'     \code{"evidence_corroborated"} (below). A fix is flagged if ANY
#'     of the following class rules fires:
#'     \itemize{
#'       \item \code{consensus}: at least 3 detectors agree
#'       \item \code{geometric_spike}: bridge AND detour
#'             (symmetric out-and-back is geometrically impossible
#'             regardless of state; the autodiff framework is blind
#'             to it, so kinematic confirmation would defeat
#'             detection)
#'       \item \code{state_anomaly}: (bridge OR detour) AND speed
#'       \item \code{kinematic_confluence}: (bridge OR detour) AND
#'             prob
#'     }
#'     Each detector contributes where its structural strengths
#'     apply, rather than being a symmetric voter. The class
#'     taxonomy used downstream (\code{error_class}) shares this
#'     structure.}
#'   \item{\code{"strict"}}{\code{(bridge OR detour) AND (prob OR
#'     speed)}.  At least one geometric AND one kinematic confirmer.
#'     Highest precision; on clean synthetic ground truth (CPF
#'     tracks) achieves zero false positives.}
#'   \item{\code{"majority"}}{Flag iff at least 2 detectors
#'     agree. Tolerates one silent detector; slightly higher recall
#'     than \code{"strict"}, slightly lower precision.}
#'   \item{\code{"speed_trusted"}}{\code{speed OR ((bridge OR detour)
#'     AND prob)}.  Trusts the speed detector alone because in its
#'     \code{"auto"} mode it is already dip-test-validated.  Catches
#'     extreme-speed transitions that strict misses, at the cost of
#'     flagging both endpoints of fast steps.}
#'   \item{\code{"competence_count"}}{Count only the detectors that are
#'     \emph{locally competent} and flag when at least \code{k} of them
#'     agree.  A detector counts at a fix only when its competence
#'     weight (supplied via \code{weight_cols}) is at least \code{gate}.
#'     With no weights every detector is fully competent everywhere and
#'     this reduces to \code{"majority"} at \code{k = 2}.  This is the
#'     mode that lets a detector's vote be silenced exactly where it is
#'     structurally blind (e.g. raw speed in the low-speed floor zone).}
#'   \item{\code{"evidence_corroborated"}}{\strong{Default.}  The
#'     evidence accumulation of \code{"weighted_evidence"} (below), but a
#'     fix is flagged only when the combined evidence is positive
#'     \emph{and} either at least two detectors corroborate (their log-LRs
#'     are positive) \emph{or} a high-specificity detector
#'     (\code{solo_cols}, default the detour detector) is \emph{saturated}
#'     -- its calibrated evidence at least 2 MADs beyond its own boundary,
#'     i.e. overwhelming within its own distribution.  This keeps the
#'     lone-strong-detector catches (e.g. a conspicuous out-and-back that
#'     only the detour detector sees) while requiring corroboration for
#'     detectors that over-react to legitimate sharp turns or
#'     behavioural-state changes (bridge, prob).  Empirically the best
#'     decision rule (benchmark 2026-06-07: highest F1 at canonical
#'     false-positive level, recall 1.0 on spoofing/jamming).  Matches the
#'     \code{\link{mt_clean_track}} cascade default.}
#'   \item{\code{"weighted_evidence"}}{Evidence accumulation. Each
#'     detector supplies a signed log-likelihood-ratio (via
#'     \code{evidence_cols}, e.g. the \code{loglr_*} columns the four
#'     primitives emit). With \code{calibrate = TRUE} (default) each
#'     detector's log-LR is first made \emph{commensurable} --
#'     \code{evidence_C * tanh(loglr / MAD(loglr))}, bounding every
#'     detector to \eqn{\pm}\code{evidence_C} so that no detector's native
#'     magnitude (e.g. the bridge's \eqn{\eta^2/2}, which can reach
#'     millions) dominates -- then summed into a combined evidence score.
#'     A fix is flagged when that score exceeds \code{evidence_threshold}
#'     (default 0, the data-driven knee between the homogeneous inlier
#'     mass and the conspicuous tail). A missing or \code{NA} contribution
#'     is an abstention (no evidence either way), \emph{not} evidence
#'     against. The combined score is returned as \code{combined_evidence}
#'     -- the single commensurable scalar behind the decision, and the
#'     substrate for an optional one-lever sensitivity control (a user
#'     threshold on that column). See \code{DESIGN_evidence_accumulation.md}.}
#'   \item{\code{"any"}}{Union of all detectors. Maximum recall,
#'     lowest precision.  Rarely the right choice for production
#'     pipelines.}
#'   \item{\code{"custom"}}{Apply a user-supplied function.  See the
#'     \code{custom} argument.}
#' }
#'
#' @param x A \code{move2} object that already has per-detector flag
#'   columns from a cascade run or from running the four primitives
#'   standalone.  Missing flag columns are treated as identically
#'   \code{FALSE} (no flags from that detector).
#' @param mode Character.  Which consensus rule to apply.  One of
#'   \code{"evidence_corroborated"} (default), \code{"class_aware"},
#'   \code{"strict"}, \code{"majority"}, \code{"speed_trusted"},
#'   \code{"competence_count"}, \code{"weighted_evidence"}, \code{"any"},
#'   or \code{"custom"}.  See \strong{Consensus modes} below.  The two
#'   evidence modes need the per-detector \code{loglr_*} columns
#'   (\code{evidence_cols}); the others need only the flag columns.
#' @param custom Function used when \code{mode = "custom"}.  Two
#'   contracts are accepted.  The legacy contract is a function of four
#'   logical vectors \code{f(by_bridge, by_prob, by_speed, by_detour)}.
#'   The general contract is a function of the evidence
#'   \code{f(flags, weights)}, where \code{flags} is a named logical
#'   matrix (one column per detector) and \code{weights} is the aligned
#'   competence matrix or \code{NULL}.  Either must return a single
#'   logical vector giving the per-fix outlier decision.  Ignored unless
#'   \code{mode = "custom"}.
#' @param detectors Named character vector mapping detector names to the
#'   flag columns that hold their decisions.  Defaults to the four
#'   standard cascade detectors.  Extend it to bring additional
#'   detectors into the decision.
#' @param weight_cols Optional named character vector mapping detector
#'   names to columns holding per-fix competence weights in \code{[0, 1]}.
#'   Used by \code{mode = "competence_count"} (and available to a custom
#'   function).  \code{NULL} (default) means every detector is fully
#'   competent everywhere, reproducing the unweighted behaviour.
#' @param evidence_cols Named character vector mapping detector names to
#'   columns holding per-fix signed log-likelihood-ratios.  Used by the
#'   \code{"evidence_corroborated"} and \code{"weighted_evidence"} modes.
#'   Defaults to the four \code{loglr_*} columns the primitives emit.  A
#'   missing column or \code{NA} entry is treated as an abstention (zero
#'   contribution); if \emph{none} of the columns are present the evidence
#'   modes error (run the detectors first, or use a flag-based mode).
#' @param evidence_threshold Numeric.  Cut for
#'   \code{mode = "weighted_evidence"}: a fix is flagged when the combined
#'   evidence exceeds this value.  Default \code{0} (the data-driven knee
#'   between the homogeneous inlier mass and the conspicuous tail; raise
#'   for a more conservative cut, lower for a more aggressive one).
#' @param calibrate Logical.  For \code{mode = "weighted_evidence"}, make
#'   the per-detector log-LRs commensurable via
#'   \code{evidence_C * tanh(loglr / MAD(loglr))} before summing.  Default
#'   \code{TRUE}; this is what keeps a high-magnitude detector from
#'   dominating the combined evidence.
#' @param evidence_C Numeric.  Saturation ceiling for the calibration:
#'   each detector contributes at most \eqn{\pm}\code{evidence_C} to the
#'   combined evidence.  Default \code{4}.
#' @param solo_cols Character vector of detector names (matching
#'   \code{evidence_cols}) allowed to flag \emph{alone} when saturated,
#'   under \code{mode = "evidence_corroborated"}.  Default \code{"detour"}
#'   -- the high-specificity geometric detector (an out-and-back is hard
#'   to produce by real movement).  Other detectors must be corroborated.
#' @param pool_by Optional grouping key for the evidence calibration -- the
#'   commensurability \code{MAD} is computed within each group, the same way
#'   the detectors pool their thresholds.  Accepts the same forms as the
#'   detector primitives: a single column name, or the nested
#'   \code{c(outer, inner)} form, looked up in \code{mt_track_data(x)} and
#'   validated by the shared resolver (a missing column errors rather than
#'   silently falling back).  Because the \code{MAD} is a distribution-scale
#'   estimate, the calibration is pooled on the \strong{outer}
#'   (distribution-source) level of a nested key; the inner union-unit role
#'   does not apply here.  \code{NULL} (default) groups per track
#'   (\code{mt_track_id}) so a multi-track object is never calibrated across
#'   tracks.  Set it to the same column(s) you pooled the detectors on.
#' @param gate Numeric in \code{[0, 1]}.  Competence threshold for
#'   \code{mode = "competence_count"}: a detector counts at a fix only
#'   when its weight is at least \code{gate}.  Default \code{0.5}.
#' @param k Integer.  Vote threshold for \code{mode =
#'   "competence_count"}: flag when at least \code{k} competent
#'   detectors agree.  Default \code{2L}.
#' @param bridge_col,prob_col,speed_col,detour_col Character.  Deprecated
#'   convenience aliases that override the corresponding entry of
#'   \code{detectors}.  Retained for backward compatibility; prefer
#'   \code{detectors}.  If a column is missing from \code{x}, that
#'   detector is treated as silent.
#'
#' @return The input \code{x} with an updated \code{is_outlier}
#'   logical column. Other per-fix columns are passed through
#'   unchanged. If \code{is_outlier} already existed on input, it
#'   is replaced. For the evidence modes (\code{"evidence_corroborated"},
#'   the default, and \code{"weighted_evidence"}) a numeric
#'   \code{combined_evidence} column is also added.
#'
#' @section Composing with custom cascades:
#'
#' If you have built your own cascade that does not use one of the
#' four standard detectors (or uses them under different column
#' names), point \code{detectors} at your column names.  Detectors with
#' no corresponding column are silently treated as \code{FALSE}, so a
#' three-detector cascade (e.g. without detour) gives the expected
#' reduced-form rules.
#'
#' For consensus rules outside the documented set, use
#' \code{mode = "custom"}.  Inside the function you can use base-R
#' logical operations, \code{rowSums(flags)} for vote counts, etc.
#'
#' @examples
#' \dontrun{
#' # Default mode (evidence_corroborated) on a cleaned cascade output
#' cleaned <- mt_clean_track(track)
#' cleaned <- mt_flag_consensus(cleaned)   # same as mt_clean_track default
#'
#' # More conservative: require a geometric AND a kinematic confirmer
#' cleaned <- mt_flag_consensus(cleaned, mode = "strict")
#'
#' # User-defined: bridge alone (geometric impossibility only)
#' cleaned <- mt_flag_consensus(cleaned, mode = "custom",
#'                              custom = function(b, p, s, d) b)
#' }
#'
#' @seealso \code{\link{mt_clean_track}} (which uses this function
#'   internally each iteration);
#'   \code{\link{mt_flag_outliers_bridge}},
#'   \code{\link{mt_flag_outliers}},
#'   \code{\link{mt_flag_speed_cap}},
#'   \code{\link{mt_flag_outliers_detour}} (the four primitives that
#'   emit the per-detector flag columns this function consumes).
#'
#' @importFrom rlang .data
#' @export
mt_flag_consensus <- function(x,
                              mode = c("evidence_corroborated", "class_aware",
                                       "strict", "majority", "speed_trusted",
                                       "competence_count", "weighted_evidence",
                                       "any", "custom"),
                              custom = NULL,
                              detectors = c(bridge = "flagged_by_bridge",
                                            prob   = "flagged_by_prob",
                                            speed  = "flagged_by_speed",
                                            detour = "flagged_by_detour"),
                              weight_cols = NULL,
                              evidence_cols = c(bridge = "loglr_bridge",
                                                prob   = "loglr_prob",
                                                speed  = "loglr_speed",
                                                detour = "loglr_detour"),
                              evidence_threshold = 0,
                              calibrate = TRUE,
                              evidence_C = 4,
                              solo_cols = "detour",
                              pool_by = NULL,
                              gate = 0.5,
                              k = 2L,
                              bridge_col = NULL,
                              prob_col   = NULL,
                              speed_col  = NULL,
                              detour_col = NULL) {
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  mode <- match.arg(mode)
  if (mode == "custom" && !is.function(custom)) {
    rlang::abort(paste0(
      "mode = \"custom\" requires `custom` to be a function with ",
      "signature `f(by_bridge, by_prob, by_speed, by_detour)` (legacy) ",
      "or `f(flags, weights)` (general) returning a logical vector."),
      class = "move2utils_mt_flag_consensus_bad_custom")
  }

  ## Backward-compat: the legacy `*_col` aliases override the
  ## corresponding entry of `detectors`.
  if (!is.null(bridge_col)) detectors["bridge"] <- bridge_col
  if (!is.null(prob_col))   detectors["prob"]   <- prob_col
  if (!is.null(speed_col))  detectors["speed"]  <- speed_col
  if (!is.null(detour_col)) detectors["detour"] <- detour_col

  n <- nrow(x)
  pull_flag <- function(col) {
    if (col %in% names(x)) {
      v <- as.logical(x[[col]])
      v[is.na(v)] <- FALSE
      v
    } else {
      rep(FALSE, n)
    }
  }
  flags <- vapply(detectors, pull_flag, logical(n))
  if (is.null(dim(flags))) {
    flags <- matrix(flags, nrow = n, dimnames = list(NULL, names(detectors)))
  }
  colnames(flags) <- names(detectors)

  weights <- NULL
  if (!is.null(weight_cols)) {
    pull_weight <- function(col) {
      if (col %in% names(x)) as.numeric(x[[col]]) else rep(1, n)
    }
    weights <- vapply(weight_cols, pull_weight, numeric(n))
    if (is.null(dim(weights))) {
      weights <- matrix(weights, nrow = n,
                        dimnames = list(NULL, names(weight_cols)))
    }
    colnames(weights) <- names(weight_cols)
  }

  evidence <- NULL
  if (!is.null(evidence_cols)) {
    pull_ev <- function(col) {
      if (col %in% names(x)) as.numeric(x[[col]]) else rep(NA_real_, n)
    }
    evidence <- vapply(evidence_cols, pull_ev, numeric(n))
    if (is.null(dim(evidence))) {
      evidence <- matrix(evidence, nrow = n,
                         dimnames = list(NULL, names(evidence_cols)))
    }
    colnames(evidence) <- names(evidence_cols)
  }

  if (mode %in% c("weighted_evidence", "evidence_corroborated")) {
    if (!any(unname(evidence_cols) %in% names(x))) {
      rlang::abort(paste0(
        "mode = \"", mode, "\" needs per-detector evidence (log-LR) ",
        "columns, but none of `evidence_cols` (",
        paste(unname(evidence_cols), collapse = ", "),
        ") are present in `x`.  Run the detectors first -- they emit the ",
        "`loglr_*` columns -- or use a flag-based mode such as ",
        "\"class_aware\"."),
        class = "move2utils_mt_flag_consensus_no_evidence")
    }
    ## Calibrate + combine on the SAME grouping the detectors used for
    ## their thresholds.  Resolve `pool_by` through the shared helper (the
    ## same one the detectors use): it looks the column up in
    ## mt_track_data(x), validates the nested c(outer, inner) form, and
    ## errors on a missing column instead of silently falling back.  The MAD
    ## is a distribution-scale estimate, so it is pooled on the OUTER
    ## (distribution-source) level; the inner union-unit role has no meaning
    ## for calibration.  NULL pool_by groups per track and never mixes tracks.
    grp <- if (!is.null(pool_by)) {
      pg <- .resolve_pool_groups(x, pool_by)
      unname(pg$outer[as.character(move2::mt_track_id(x))])
    } else {
      as.character(move2::mt_track_id(x))
    }
    Cm <- .combine_evidence(evidence, calibrate = calibrate, C = evidence_C,
                            group = grp, components = TRUE)
    E <- rowSums(Cm)
    ## Expose the combined evidence -- the single commensurable scalar
    ## behind the decision, and the substrate for the one-lever
    ## sensitivity UI (a user threshold on this column).
    x$combined_evidence <- E
    if (mode == "weighted_evidence") {
      x$is_outlier <- E > evidence_threshold
    } else {
      ## evidence_corroborated (default): positive net evidence AND either
      ## >=2 detectors corroborate, OR a high-specificity solo detector
      ## (`solo_cols`, default detour) is saturated -- its calibrated
      ## evidence >= 2 MADs beyond its own boundary (overwhelming in its
      ## own distribution).  Bridge/prob over-react to legitimate sharp
      ## turns / behavioural-state changes, so a lone one must be
      ## corroborated; an out-and-back (detour) is unambiguous and may
      ## carry alone.  Matches the mt_clean_track cascade rule.
      pos <- rowSums(evidence > 0, na.rm = TRUE)
      solo_idx <- which(colnames(Cm) %in% solo_cols)
      sat <- if (length(solo_idx)) {
        apply(Cm[, solo_idx, drop = FALSE], 1L, max) > evidence_C * tanh(2)
      } else rep(FALSE, nrow(Cm))
      x$is_outlier <- (E > evidence_threshold) & (pos >= 2L | sat)
    }
  } else {
    x$is_outlier <- .consensus_decide(flags, weights, mode = mode,
                                      custom = custom, gate = gate, k = k,
                                      evidence = evidence,
                                      evidence_threshold = evidence_threshold)
  }
  x
}


## Internal core: the actual logic, operating on a named logical flag
## matrix (one column per detector) plus an optional aligned numeric
## competence-weight matrix.  Exposed separately so mt_clean_track can
## call it from inside its iteration loop without re-materialising
## columns on the move2.
##
## Behaviour invariant: with the four standard detector columns and
## NULL weights, every legacy mode is byte-identical to the pre-refactor
## four-vector implementation (see tests/testthat/test-consensus-contract.R).
##
## Combine per-detector evidence columns (signed log-LRs) into a single
## posterior score.  With calibrate = TRUE each column is first made
## commensurable via C * tanh(loglr / MAD) -- bounding every detector to
## +/- C so that no detector's native magnitude (e.g. the bridge's
## eta^2/2, which can reach millions) dominates the sum -- then summed.
## NA (abstention) contributes 0.  The MAD scale is taken per detector
## over its finite values on this object (self-calibrating, unsupervised).
## @keywords internal
## When `components = TRUE`, return the per-detector calibrated matrix
## (one column per detector) instead of its row sums -- used by the
## evidence_corroborated rule, which needs both the sum and the
## individual saturated columns.
.combine_evidence <- function(evidence, calibrate = FALSE, C = 4,
                              group = NULL, components = FALSE) {
  n <- if (is.null(evidence)) 0L else nrow(evidence)
  if (is.null(evidence) || ncol(evidence) == 0L) {
    return(if (components) matrix(0, n, 0L) else rep(0, n))
  }
  if (!calibrate) {
    ev <- evidence
    ev[is.na(ev)] <- 0
    return(if (components) ev else rowSums(ev))
  }
  cal1 <- function(v) {
    sc <- stats::mad(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(sc) || sc == 0) sc <- 1
    o <- C * tanh(v / sc)
    o[is.na(o)] <- 0
    o
  }
  cols <- colnames(evidence)
  if (is.null(cols)) cols <- paste0("V", seq_len(ncol(evidence)))
  if (is.null(group) || length(unique(group)) <= 1L) {
    ev <- vapply(seq_len(ncol(evidence)),
                 function(j) cal1(evidence[, j]), numeric(n))
    if (is.null(dim(ev))) ev <- matrix(ev, nrow = n)
    colnames(ev) <- cols
    return(if (components) ev else rowSums(ev))
  }
  ## group-aware: calibrate each column within its group so the MAD scale
  ## is pooled on the same key the detectors used for their thresholds.
  ## dplyr group_by (preserves row order under mutate).
  d <- as.data.frame(evidence)
  names(d) <- cols
  d$.grp <- as.character(group)
  d <- dplyr::group_by(d, dplyr::across(dplyr::all_of(".grp")))
  d <- dplyr::mutate(d, dplyr::across(dplyr::all_of(cols), cal1))
  d <- dplyr::ungroup(d)
  ev <- as.matrix(d[cols])
  if (components) ev else rowSums(ev)
}

## @keywords internal
.consensus_decide <- function(flags, weights = NULL,
                              mode = "class_aware", custom = NULL,
                              gate = 0.5, k = 2L,
                              evidence = NULL, evidence_threshold = 0,
                              calibrate = FALSE, evidence_C = 4) {
  n  <- nrow(flags)
  cn <- colnames(flags)
  col <- function(nm) if (!is.null(cn) && nm %in% cn) flags[, nm] else rep(FALSE, n)
  bridge <- col("bridge"); prob   <- col("prob")
  speed  <- col("speed");  detour <- col("detour")
  geo    <- bridge | detour
  votes  <- if (ncol(flags) > 0L) rowSums(flags) else rep(0L, n)

  switch(mode,
    class_aware = {
      ## Union of class-specific rules.  Each rule encodes where each
      ## detector is in scope; single-detector fires never flag.
      fire_consensus  <- votes >= 3L
      fire_geom_spike <- bridge & detour
      fire_state_anom <- geo & speed
      fire_kin_conf   <- geo & prob
      fire_consensus | fire_geom_spike | fire_state_anom | fire_kin_conf
    },
    strict        = geo & (prob | speed),
    majority      = votes >= 2L,
    speed_trusted = speed | (geo & prob),
    any           = votes >= 1L,
    competence_count = {
      ## Count only detectors that are locally competent (weight >= gate)
      ## and require at least k of them to agree.  NULL weights => every
      ## detector competent everywhere => reduces to majority at k = 2.
      walign <- matrix(1, nrow = n, ncol = ncol(flags),
                       dimnames = list(NULL, cn))
      if (!is.null(weights)) {
        common <- intersect(cn, colnames(weights))
        if (length(common)) walign[, common] <- weights[, common]
      }
      rowSums(flags & (walign >= gate)) >= k
    },
    weighted_evidence = {
      ## Evidence accumulation: sum signed per-detector log-LRs into a
      ## posterior log-odds and threshold it.  Each column is one
      ## detector's competence-scaled log-likelihood-ratio of outlier vs
      ## not.  NA (or a missing detector) means abstain = LR 1 = 0
      ## contribution -- NOT evidence against.  The three LR kinds
      ## (two-sided density-ratio, one-sided hard-constraint, folded
      ## correlated streams) are encoded upstream in how each column is
      ## filled; the accumulator only sums.
      if (is.null(evidence)) {
        rlang::abort(paste0(
          "mode = \"weighted_evidence\" requires an `evidence` matrix of ",
          "per-detector log-likelihood-ratios."),
          class = "move2utils_mt_flag_consensus_missing_evidence")
      }
      .combine_evidence(evidence, calibrate = calibrate, C = evidence_C) >
        evidence_threshold
    },
    custom = {
      ## Two accepted contracts.  A four-argument function is the legacy
      ## f(by_bridge, by_prob, by_speed, by_detour); anything else is the
      ## general f(flags, weights).
      if (length(formals(custom)) == 4L) {
        out <- custom(bridge, prob, speed, detour)
      } else {
        out <- custom(flags, weights)
      }
      if (!is.logical(out) || length(out) != n) {
        rlang::abort(paste0(
          "`custom` must return a logical vector of the same length ",
          "as its inputs."),
          class = "move2utils_mt_flag_consensus_bad_custom_return")
      }
      out[is.na(out)] <- FALSE
      out
    }
  )
}


## Thin backward-compatible wrapper over .consensus_decide() preserving
## the original four-vector signature.  mt_clean_track() calls this from
## its iteration loop, so its signature must not change.
##
## @keywords internal
.consensus_logical <- function(by_bridge, by_prob, by_speed, by_detour,
                                mode = "class_aware", custom = NULL) {
  flags <- cbind(bridge = as.logical(by_bridge),
                 prob   = as.logical(by_prob),
                 speed  = as.logical(by_speed),
                 detour = as.logical(by_detour))
  .consensus_decide(flags, weights = NULL, mode = mode, custom = custom)
}

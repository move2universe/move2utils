#' @keywords internal
#' @aliases move2utils-package
#'
#' @section Outlier detection:
#' The scientific core is a set of detectors plus a one-call orchestrator.
#' All accept `move2` objects in any CRS (longitude/latitude is
#' auto-projected to a local AEQD internally; output is returned in the
#' caller's CRS) and handle single- or multi-track input per individual.
#' \itemize{
#'   \item [mt_clean_track()] -- one-call entry point: orchestrates the four
#'     detectors with an evidence-accumulation decision rule and topological
#'     block expansion, iterating to convergence.
#'   \item [mt_flag_outliers()] -- movement-metric probability detector.
#'   \item [mt_flag_outliers_bridge()] -- Brownian-bridge geometric residual.
#'   \item [mt_flag_outliers_detour()] -- path-vs-displacement (time-insensitive).
#'   \item [mt_flag_speed_cap()] -- step-level speed cap; companion
#'     [mt_suggest_speed_cap()].
#'   \item [mt_sequential_outliers()], [mt_combined_outliers()],
#'     [mt_peel_speed()] -- sequential / voting / iterative-peel variants.
#'   \item [mt_flag_outliers_dbgb()] / [mt_flag_outliers_dbbmm()] -- state-aware
#'     bridge detectors (standalone diagnostic).
#'   \item [mt_flag_consensus()] -- turn several detectors' flags into one
#'     decision; [mt_persistence_score()], [mt_diagnose_clean_track()],
#'     [mt_diagnose_flags()] -- annotation / post-run diagnostics.
#'   \item [v_phys_estimate()] -- allometric maximum speed from mass + mode.
#' }
#'
#' @section Utilisation distributions and motion variance:
#' [mt_dbbmm_variance()] / [mt_dbgb_variance()] (isotropic / directional
#' dynamic Brownian-bridge variance), [mt_dbbmm_ud()] / [mt_dbgb_ud()]
#' (utilisation distributions), [mt_motion_variance()],
#' [mt_suggest_dbbmm_window()], and the raster utilities [ud_volume()],
#' [ud_outer_probability()], [emd()]. These accept any CRS and return the UD
#' in the caller's CRS.
#'
#' @section Track utilities:
#' [mt_corridor()], [mt_thin_time()], [mt_thin_distance()],
#' [mt_mask_segments()], [mt_filter_gps_quality()].
#'
#' @section Threshold-type vocabulary:
#' The data-driven detectors expose a `threshold_type` selecting how the
#' flag boundary is found from the data:
#' \itemize{
#'   \item `"gap"` -- broken-stick gap on the sorted log-scores.
#'   \item `"entropy"` -- density-valley (entropy) cut; safe on clean,
#'     unimodal data.
#'   \item `"significance"` / `"percentile"` -- parametric / quantile cuts
#'     (probability detector).
#'   \item `"auto"` -- self-determining cut (detour, speed cap).
#'   \item `"fixed"` / `"hard"` -- a user-supplied literal threshold / cap.
#' }
#'
#' @section Consensus vocabulary:
#' [mt_flag_consensus()] (and `mt_clean_track(consensus = )`) turn the
#' detectors' flags into a decision. The default is
#' `"evidence_corroborated"`: each detector emits a signed log-likelihood
#' ratio, these are calibrated commensurably into a `combined_evidence`
#' score, and a fix is flagged when the evidence is positive AND either two
#' detectors corroborate or the high-specificity detour detector is
#' saturated. Other modes: `"weighted_evidence"` (the smooth evidence score
#' thresholded directly), the Boolean `"class_aware"` (the prior default),
#' `"strict"`, `"majority"`, `"any"`, `"competence_count"`, and `"custom"`.
#'
#' @section Conventions:
#' Functions are prefixed `mt_` when they take `move2` objects (raster /
#' scalar utilities are not). Preconditions are signalled with stable classed
#' conditions (`move2utils_<context>_<reason>`); empty or non-finite geometry
#' is rejected rather than silently dropped. See `vignette(package =
#' "move2utils")` for worked examples.
"_PACKAGE"

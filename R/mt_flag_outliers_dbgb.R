## State-aware bridge-residual outlier detection via dBBMM-derived
## motion variance.
##
## The legacy bridge primitive (mt_flag_outliers_bridge) normalises
## residuals by sqrt(width^2) only -- the time-only Brownian-bridge
## SD when sigma^2 is set to 1. That makes it leverage-immune (no
## variance estimated from the data) but also state-blind: on a track
## with multiple behavioural states the legacy eta distribution
## becomes multimodal, and the threshold detector lands in the valley
## between activity modes rather than above the rightmost mode.
##
## This primitive uses dBBMM/dBGB-derived sigma^2(t) to normalise
## residuals into a state-aware Z-score:
##
##   Z(i) = residual(i) / sqrt(sigma^2(t_i) * width(i)^2 + 2 * location_error^2)
##
## Z follows Rayleigh(1) under the Brownian-bridge null, so the
## threshold is principled: P(Z > z) = exp(-z^2 / 2).
##
## Trade-off: state-aware but exposed to the leverage problem (extreme
## outliers inflate sigma^2 estimation and absorb themselves into
## high-sigma^2 segments).  Mitigation: run an upstream speed peel
## (Stage 1 of the cascade) before calling this primitive.  When
## pre_peel_v_max is supplied the primitive does the peel internally;
## otherwise it is the caller's responsibility.


#' Flag outliers via state-aware bridge residuals (dBGB envelope)
#'
#' State-aware bridge-residual primitive.  Estimates per-axis
#' (parallel / orthogonal) motion variance via
#' \code{\link{mt_dbgb_variance}} and flags fixes by an envelope rule
#' across three Z channels: along-axis (\eqn{Z_{\parallel}}),
#' across-axis (\eqn{Z_{\perp}}), and the joint
#' \eqn{Z_{\chi^2_2} = \sqrt{Z_{\parallel}^2 + Z_{\perp}^2}}.  Under
#' the constraint \eqn{\sigma^2_{\parallel} = \sigma^2_{\perp}} (i.e.
#' isotropic motion variance) the dBGB model reduces to dBBMM and
#' \eqn{Z_{\chi^2_2}} reduces to the dBBMM scalar-Z; that submodel is
#' available via the \code{variance = "dbbmm"} option here, or
#' equivalently via \code{\link{mt_flag_outliers_dbbmm}} as a
#' speed-saving wrapper.  The leverage-immune static bridge primitive
#' is \code{\link{mt_flag_outliers_bridge}}.
#'
#' @details
#' \strong{Three Z channels.}  Under \code{variance = "dbgb"}, the
#' bridge residual at each fix is decomposed onto the local travel
#' axis (\eqn{r_{\parallel}}) and its perpendicular
#' (\eqn{r_{\perp}}).  Each component is normalised by its own
#' axis-specific bridge variance:
#'
#'   \deqn{Z_{\parallel,i} = r_{\parallel,i} /
#'         \sqrt{\sigma^2_{\parallel}(t_i) \cdot w_i^2 +
#'                2 \cdot \mathrm{location\_error}^2}}
#'
#' and analogously for \eqn{Z_{\perp,i}}.  The chisq channel is
#' \eqn{Z_{\chi^2_2,i} = \sqrt{Z_{\parallel,i}^2 + Z_{\perp,i}^2}},
#' which under the diagonal-Sigma bridge null is Rayleigh(1)-
#' distributed -- the same null distribution as the dBBMM scalar Z.
#'
#' \strong{Envelope flag rule.}  A fix is flagged when ANY of the
#' three Z channels exceeds its own threshold.  This catches three
#' distinct outlier signatures:
#' \itemize{
#'   \item \eqn{Z_{\parallel}}-only: a step-acceleration anomaly along
#'     the local travel direction (overshoot / undershoot in the
#'     direction of motion).
#'   \item \eqn{Z_{\perp}}-only: a sideways anomaly perpendicular to
#'     travel (typical of colony-return spikes in migrating animals
#'     or azimuthal GPS jitter).
#'   \item \eqn{Z_{\chi^2_2}}-only: a combined-axis anomaly where
#'     neither component exceeds its per-axis threshold but their
#'     quadrature sum does (an isotropic spike with no axis
#'     preference).
#' }
#'
#' \strong{Threshold methods.}  Three calibrations are available via
#' \code{z_threshold_method}:
#' \itemize{
#'   \item \code{"bonferroni"} (default).  Per-channel Bonferroni-Z
#'     thresholds with joint FWER \eqn{\le 0.05}: \eqn{\alpha/2} to
#'     the chisq channel (Rayleigh-tail), \eqn{\alpha/4} to each per-
#'     axis channel (\eqn{|N(0,1)|}-tail).  Conservative and FWER-
#'     correct.
#'   \item \code{"bh_fdr"}.  Same alpha allocation but using
#'     Benjamini-Hochberg FDR per channel.  More sensitive at the
#'     same nominal alpha; appropriate when missing a real outlier
#'     costs more than flagging a marginal one.
#'   \item \code{"gap"}.  Per-channel data-driven break via the
#'     package's broken-stick + tail-decay inflection helper.  The
#'     same self-thresholding family used elsewhere in the package
#'     (e.g. \code{\link{mt_flag_outliers_bridge}}, the
#'     \code{residual_max} gate below).  Note: the gap detector was
#'     calibrated for log-probability scores -- on raw Z (whose null
#'     is half-normal / Rayleigh with thin natural tails) it tends to
#'     find a "break" at the natural decay edge of the bulk and
#'     over-flags clean tracks.  Use only when the parametric
#'     calibration is too conservative for your application.
#' }
#'
#' \strong{Diagnostic labelling.}  Every flagged fix carries a
#' \code{bridge_z_class} label identifying which channel(s) fired
#' and (for chisq-only flags) which axis dominated.  This is the
#' principal value-add over the dBBMM-only primitive: the flag set
#' is typically the same on real tracks (the per-axis envelope
#' tests rarely fire under proper FWER control), but the directional
#' attribution lets users explain \emph{why} each fix was caught.
#'
#' @section The leverage caveat:
#' The dBBMM variance estimator is sensitive to outliers (large
#' residuals inflate the windowed sigma^2 estimate, after which the
#' offending fix's Z score is bounded -- the leverage problem).  The
#' legacy bridge avoids this by not estimating variance.  If your
#' track contains physiologically-impossible-speed events, run
#' \code{\link{mt_peel_speed}} first OR pass \code{pre_peel_v_max}
#' (which runs the peel internally before sigma^2 estimation).
#' \code{\link{mt_clean_track}} composes the two stages by default.
#'
#' @param x A \code{move2} object.  Single- or multi-track.  Auto-
#'   projected to a per-track local AEQD if input is lon/lat.
#' @param location_error Per-fix GPS noise floor in metres (1-sigma),
#'   accepted in the same forms as the other primitives: \code{NULL}
#'   (default -- no anchor noise floor), a numeric scalar, a per-fix
#'   numeric vector of length \code{nrow(x)}, a column name, or the string
#'   \code{"auto"} (reads Movebank quality columns).  Resolved to a numeric
#'   vector internally.  Enters the denominator as
#'   \code{2 * location_error^2}, the variance contribution from
#'   imperfect knowledge of both bridge anchors; with the default, that
#'   term is zero and the denominator reduces to
#'   \code{sqrt(sigma^2(t_i) * w_i^2)}.  For e-obs / Argos-A class
#'   hardware, \code{25} is a defensible scalar; for other devices,
#'   supply a calibration estimate.  See the package vocabulary note:
#'   \code{location_error} is a \emph{user-supplied} calibration prior,
#'   not a Movebank-exposed device attribute -- those are handled
#'   separately as quality columns.
#' @param z_threshold Numeric scalar or \code{NULL}.  When supplied,
#'   overrides only the chisq-channel threshold; per-axis thresholds
#'   stay parametric.  Incompatible with
#'   \code{z_threshold_method = "gap"} (which is data-driven).
#' @param z_threshold_method Character, one of \code{"bonferroni"}
#'   (default), \code{"bh_fdr"}, or \code{"gap"}.  Selects how the
#'   per-channel thresholds are derived; see Details for the alpha
#'   allocation across the three Z channels.
#' @param residual_max Numeric scalar, \code{Inf}, or \code{NULL}.
#'   The absolute-residual axis: a fix is flagged if its bridge
#'   residual exceeds this magnitude regardless of \eqn{Z}.  Catches
#'   geometrically extreme residuals that dBBMM-Z would absorb into
#'   local sigma^2.  Default \code{NULL} runs a "last-satisfactory-
#'   break" gap detector on \code{log(residual)}: it walks gaps from
#'   the extreme tail down toward the bulk and picks the *last* gap
#'   that exceeds the broken-stick null.  This places the cap at the
#'   bulk-to-tail boundary, ensuring all residuals above the bulk are
#'   caught even when the outlier set itself is multi-scale (e.g.
#'   spoof clusters at multiple distance scales).  Track-specific:
#'   stork bulk ~100 m gets cap ~ 1 km; swift bulk in km range gets
#'   cap scaled accordingly.  \code{Inf} disables the axis.
#' @param window_size Integer or \code{NULL}.  dBBMM sliding-window
#'   size.  \code{NULL} (default) auto-selects from the data's motion
#'   autocorrelation: window = \code{c} * decorrelation lag of
#'   \eqn{v^2}, with \code{c = 4}.  Capped at \code{max(11, n / 8)}.
#' @param margin Integer or \code{NULL}.  dBBMM margin.  \code{NULL}
#'   (default) sets \code{margin = (window_size - 1) / 4}, rounded
#'   to the nearest odd integer.
#' @param z_gap_threshold Numeric multiplier controlling strictness of
#'   the broken-stick gap detector when
#'   \code{z_threshold_method = "gap"}.  \code{NULL} (default) defers
#'   to \code{.z_gap_threshold}'s leaf formal (3).  Higher = more
#'   conservative (fewer flags).  Plausible range 2--5.
#' @param residual_entropy_threshold,residual_gap_threshold Numeric
#'   scalars or \code{NULL}.  Control the entropy / gap detectors
#'   inside the \code{residual_max = NULL} auto path -- the gate
#'   on \code{|residual|} in metres that catches geometrically
#'   extreme residuals.  \code{NULL} (default) defers to the leaf
#'   formals.  See audits/2026-05-25-parameter-propagation/findings.md
#'   §1.2.
#' @param residual_dip_alpha Numeric in (0, 1).  Significance level
#'   for the Hartigan dip test that validates the broken-stick
#'   fallback inside the \code{residual_max} auto path.  Default
#'   \code{0.05} (Fisher convention).
#' @param variance Character, one of \code{"dbgb"} (default,
#'   directional) or \code{"dbbmm"} (isotropic special case, faster
#'   but discards the directional decomposition).  Under
#'   \code{"dbbmm"} only the chisq channel exists, and the envelope
#'   reduces to the legacy scalar Rayleigh-Z test.
#' @param pre_peel_v_max Numeric scalar or \code{NULL}.  If supplied,
#'   run \code{\link{mt_peel_speed}} at this cap before sigma^2
#'   estimation; this protects against the leverage problem when
#'   physiologically-impossible fixes are present in the data.
#'   \code{NULL} (default) skips the peel -- standalone use of this
#'   primitive.
#' @param plot Logical.  Default TRUE; diagnostic plot.
#' @param remove Logical.  If \code{TRUE}, drop flagged rows from the
#'   returned object.  Default \code{FALSE}.
#' @param silent Logical.  Default \code{FALSE}.
#'
#' @return The input object with added columns:
#'   \describe{
#'     \item{\code{bridge_residual}}{Euclidean residual from bridge
#'           mean, in metres.}
#'     \item{\code{bridge_width}}{Bridge width (sqrt-time scale).}
#'     \item{\code{sigma2_motion_para},
#'           \code{sigma2_motion_orth}}{Per-axis local
#'           \eqn{\sigma^2} from dBGB in m^2/min.  Only emitted under
#'           \code{variance = "dbgb"}.}
#'     \item{\code{sigma2_motion}}{Isotropic local
#'           \eqn{\sigma^2} from dBBMM in m^2/min.  Only emitted under
#'           \code{variance = "dbbmm"}.}
#'     \item{\code{bridge_z}}{The chisq Z used as the principal
#'           flagging statistic (= \code{bridge_z_chisq}).  Kept for
#'           back-compat with downstream code that reads a single
#'           \code{bridge_z} column.}
#'     \item{\code{bridge_z_chisq}}{The combined-axis Z,
#'           \eqn{\sqrt{Z_{\parallel}^2 + Z_{\perp}^2}}.  Rayleigh(1)
#'           under the bridge null.  Equal to \code{bridge_z} under
#'           \code{variance = "dbbmm"}.}
#'     \item{\code{bridge_z_para}, \code{bridge_z_orth}}{Per-axis
#'           half-normal Z-scores from the dBGB decomposition.  Only
#'           emitted under \code{variance = "dbgb"}.}
#'     \item{\code{bridge_z_class}}{Diagnostic label classifying which
#'           channel(s) fired on flagged fixes (\code{"para"},
#'           \code{"orth"}, \code{"isotropic"},
#'           \code{"isotropic_para_dom"} / \code{"isotropic_orth_dom"}
#'           sub-classes for chisq-only flags by axis dominance,
#'           multi-channel combinations, \code{"residual"} for the
#'           geometric-impossibility cap, or \code{"none"}).}
#'     \item{\code{bridge_iteration}}{Integer iteration at which the
#'           fix was flagged (always 0 / 1 here -- this primitive is
#'           single-pass; provided for API symmetry).}
#'     \item{\code{is_outlier}}{Logical.  TRUE where flagged.}
#'   }
#'
#'   Plus attributes \code{z_thresholds} (list with per-channel
#'   thresholds), \code{z_threshold} (the chisq threshold, for back-
#'   compat), \code{z_threshold_method}, \code{residual_max},
#'   \code{window_size}, \code{location_error}.
#'
#' @examples
#' \dontrun{
#' library(move2)
#' x <- mt_read(system.file("extdata/synthetic_tracks.csv.gz",
#'                            package = "move2utils"))
#' x <- filter_track_data(x, .track_id = "CPF_A")
#'
#' ## standalone (single-state track or pre-cleaned)
#' res <- mt_flag_outliers_dbgb(x)
#'
#' ## cascade composition: peel impossible-speed fixes first
#' res <- mt_flag_outliers_dbgb(x, pre_peel_v_max = 50)
#' }
#'
#' @seealso \code{\link{mt_flag_outliers_dbbmm}} (isotropic
#'   submodel; speed-saving wrapper around this function with
#'   \code{variance = "dbbmm"});
#'   \code{\link{mt_flag_outliers_bridge}} (leverage-immune static
#'   bridge); \code{\link{mt_peel_speed}} (Stage 1 of the cascade);
#'   \code{\link{mt_dbgb_variance}}; \code{\link{mt_clean_track}}.
#'
#' @importFrom move2 mt_time mt_track_id mt_aeqd_crs mt_n_tracks
#' @importFrom sf st_coordinates st_is_longlat st_transform
#' @importFrom stats acf
#' @importFrom graphics par plot abline points lines legend
#' @export
mt_flag_outliers_dbgb <- function(x,
                                    location_error    = NULL,
                                    z_threshold       = NULL,
                                    z_threshold_method = c("bonferroni",
                                                            "bh_fdr",
                                                            "gap"),
                                    z_gap_threshold    = NULL,
                                    residual_max      = NULL,
                                    residual_entropy_threshold = NULL,
                                    residual_gap_threshold     = NULL,
                                    residual_dip_alpha         = 0.05,
                                    window_size       = NULL,
                                    margin            = NULL,
                                    variance          = c("dbgb", "dbbmm"),
                                    pre_peel_v_max    = NULL,
                                    plot              = TRUE,
                                    remove            = FALSE,
                                    silent            = FALSE) {
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  variance           <- match.arg(variance)
  z_threshold_method <- match.arg(z_threshold_method)
  say <- .say(silent)

  ## ---- multi-individual dispatch ----
  ids <- move2::mt_track_id(x)
  unique_ids <- unique(ids)
  if (length(unique_ids) > 1L) {
    say("Processing ", length(unique_ids), " individuals separately...")
    parts <- lapply(unique_ids, function(id) {
      xi <- x[ids == id, ]
      say(sprintf("--- %s (%d locations) ---", id, nrow(xi)))
      mt_flag_outliers_dbgb(xi,
        location_error = location_error, z_threshold = z_threshold,
        z_threshold_method = z_threshold_method,
        z_gap_threshold = z_gap_threshold,
        residual_max = residual_max,
        residual_entropy_threshold = residual_entropy_threshold,
        residual_gap_threshold     = residual_gap_threshold,
        residual_dip_alpha         = residual_dip_alpha,
        window_size = window_size, margin = margin,
        variance = variance,
        pre_peel_v_max = pre_peel_v_max,
        plot = FALSE, remove = FALSE, silent = silent)
    })
    out <- do.call(rbind, parts)
    if (remove) out <- out[!out$is_outlier, ]
    return(out)
  }

  ## ---- auto-project lon/lat ----
  orig_was_longlat <- isTRUE(sf::st_is_longlat(x))
  if (orig_was_longlat) {
    say("Input is in longitude/latitude. Auto-projecting to local AEQD; ",
        "output is returned in the original CRS.")
    orig_x <- x
    aeqd   <- move2::mt_aeqd_crs(x, center = "center", units = "m")
    x      <- sf::st_transform(x, aeqd)
  }

  n_total <- nrow(x)
  if (n_total < 10L) {
    rlang::warn("Too few locations (<10). Returning without flags.",
                class = "move2utils_mt_flag_outliers_dbgb_too_few_locations")
    x$bridge_residual   <- rep(NA_real_, n_total)
    x$bridge_width      <- rep(NA_real_, n_total)
    x$sigma2_motion     <- rep(NA_real_, n_total)
    x$bridge_z          <- rep(NA_real_, n_total)
    x$bridge_z_chisq    <- rep(NA_real_, n_total)
    x$bridge_z_class    <- rep(NA_character_, n_total)
    if (variance == "dbgb") {
      x$bridge_z_para <- rep(NA_real_, n_total)
      x$bridge_z_orth <- rep(NA_real_, n_total)
    }
    x$bridge_iteration  <- rep(NA_integer_, n_total)
    x$is_outlier        <- rep(FALSE, n_total)
    if (orig_was_longlat) {
      cols <- c("bridge_residual","bridge_width","sigma2_motion",
                "bridge_z","bridge_z_chisq","bridge_z_class",
                "bridge_iteration","is_outlier")
      if (variance == "dbgb") cols <- c(cols, "bridge_z_para", "bridge_z_orth")
      for (col in cols) orig_x[[col]] <- x[[col]]
      x <- orig_x
    }
    if (remove) x <- x[!x$is_outlier, ]
    return(x)
  }

  ## ---- Stage 1: optional pre-peel for leverage protection ----
  is_peeled <- rep(FALSE, n_total)
  if (!is.null(pre_peel_v_max)) {
    if (!is.numeric(pre_peel_v_max) || length(pre_peel_v_max) != 1L ||
        is.na(pre_peel_v_max) || pre_peel_v_max <= 0)
      rlang::abort("`pre_peel_v_max` must be a positive scalar (m/s) or NULL.",
                   class = "move2utils_mt_flag_outliers_dbgb_bad_pre_peel_v_max")
    pp <- suppressMessages(
      mt_peel_speed(x, v_max = pre_peel_v_max, remove = FALSE))
    is_peeled <- pp$is_outlier
    say(sprintf("Stage 1 (pre-peel) at v_max = %g m/s: %d / %d removed.",
                pre_peel_v_max, sum(is_peeled), n_total))
  }

  ## ---- Stage 2: sigma^2 estimation on survivors ----
  surv_idx <- which(!is_peeled)
  surv_x   <- x[surv_idx, ]
  ns       <- nrow(surv_x)

  ## Resolve location_error ONCE to a clean per-fix numeric vector aligned to
  ## x (accepts NULL / scalar / per-fix vector / column-name / "auto", like
  ## the other primitives), then slice to the survivors.  This resolves
  ## "auto"/column inputs to numbers *before* they reach the 2*location_error^2
  ## Z-denominator (which would error on a character) and aligns a per-fix
  ## vector to surv_x rather than the full track.  NA error -> 0 (no inflation).
  le_full <- .resolve_location_error(location_error, x, n_total)
  if (is.null(le_full)) le_full <- rep(0, n_total)
  le_full[is.na(le_full)] <- 0
  le_surv <- le_full[surv_idx]

  ## resolve window_size + margin
  if (is.null(window_size)) {
    window_size <- .auto_window_size_acf(surv_x)
    say(sprintf("Auto window_size from ACF(v^2): %d", window_size))
  }
  ## hard caps on window_size
  window_size <- max(11L, as.integer(window_size))
  if (window_size %% 2L == 0L) window_size <- window_size - 1L
  window_size <- min(window_size, max(11L, as.integer(ns / 8L)))
  if (window_size %% 2L == 0L) window_size <- window_size - 1L
  if (is.null(margin)) {
    margin <- max(5L, as.integer((window_size - 1L) / 4L))
    if (margin %% 2L == 0L) margin <- margin - 1L
  }

  if (ns < 2L * window_size) {
    rlang::warn(
      sprintf("Only %d survivors after Stage 1; below 2 * window_size (%d). Returning without dBBMM flags.",
              ns, 2L * window_size),
      class = "move2utils_mt_flag_outliers_dbgb_too_few_survivors")
    x$bridge_residual   <- rep(NA_real_, n_total)
    x$bridge_width      <- rep(NA_real_, n_total)
    x$sigma2_motion     <- rep(NA_real_, n_total)
    x$bridge_z          <- rep(NA_real_, n_total)
    x$bridge_z_chisq    <- rep(NA_real_, n_total)
    x$bridge_z_class    <- rep(NA_character_, n_total)
    if (variance == "dbgb") {
      x$bridge_z_para <- rep(NA_real_, n_total)
      x$bridge_z_orth <- rep(NA_real_, n_total)
    }
    x$bridge_iteration  <- rep(NA_integer_, n_total)
    x$is_outlier        <- is_peeled
    if (orig_was_longlat) {
      cols <- c("bridge_residual","bridge_width","sigma2_motion",
                "bridge_z","bridge_z_chisq","bridge_z_class",
                "bridge_iteration","is_outlier")
      if (variance == "dbgb") cols <- c(cols, "bridge_z_para", "bridge_z_orth")
      for (col in cols) orig_x[[col]] <- x[[col]]
      x <- orig_x
    }
    if (remove) x <- x[!x$is_outlier, ]
    return(x)
  }

  variance_fun <- if (variance == "dbgb")
                    move2utils::mt_dbgb_variance
                  else
                    move2utils::mt_dbbmm_variance
  vo <- variance_fun(surv_x, location_error = le_surv,
                      window_size = window_size, margin = margin)
  mv <- mt_motion_variance(vo)
  say(sprintf("Stage 2: %s sigma^2 estimated (window=%d, margin=%d, n=%d).",
              variance, window_size, margin, ns))

  ## ---- Stage 3: residual + Z + threshold ----
  ## Two physical Z channels under variance = "dbgb":
  ##   - Z_para, Z_orth: per-axis residual / per-axis sigma_eff.
  ##     Each is |Normal(0,1)|-distributed under the bridge null.
  ##   - Z_chisq = sqrt(Z_para^2 + Z_orth^2): the 2D Mahalanobis
  ##     magnitude under diagonal Sigma.  Rayleigh(1) under the null.
  ##     Reduces to the dBBMM scalar-Z when sigma_para = sigma_orth.
  ## The three signals carve out three rejection regions in the
  ## (Z_para, Z_orth) plane: two axis-aligned strips and an outside-
  ## circle.  The envelope rule fires if any of the three exceeds
  ## its threshold and labels each flagged fix with which region(s)
  ## it landed in (bridge_z_class).  Under variance = "dbbmm" the
  ## per-axis decomposition is unavailable; the function falls back
  ## to the scalar Rayleigh-Z path with no envelope.
  if (variance == "dbgb") {
    bridge <- .compute_bridge_z_dbgb(surv_x, mv, le_surv)
    bridge$Z_chisq <- sqrt(bridge$Z_para ^ 2 + bridge$Z_orth ^ 2)
    bridge$Z <- bridge$Z_chisq
  } else {
    bridge <- .compute_bridge_z(surv_x, mv, le_surv)
    bridge$Z_para  <- rep(NA_real_, length(bridge$Z))
    bridge$Z_orth  <- rep(NA_real_, length(bridge$Z))
    bridge$Z_chisq <- bridge$Z
  }

  ## Resolve per-channel thresholds.  All three modes produce one
  ## threshold per channel and the envelope flag rule fires if any
  ## channel exceeds its own threshold.
  ##
  ## - "bonferroni" (default): per-channel Bonferroni-Z thresholds
  ##   under joint FWER <= alpha = 0.05.  Allocation: alpha/2 to the
  ##   chisq channel (Rayleigh null on Z_chisq, one test family of
  ##   ns fixes) and alpha/2 split equally between the two per-axis
  ##   channels (|N(0,1)| null, alpha/4 per axis over ns fixes each).
  ##   Joint FWER <= alpha by Bonferroni union.  Conservative; the
  ##   tightest calibration under the bridge null and the right
  ##   default when controlling false-positive rate matters.
  ## - "bh_fdr": same alpha allocation but using BH-FDR per channel
  ##   (Rayleigh-tail BH on Z_chisq, |N(0,1)|-tail BH on each axis).
  ##   More sensitive at controlled FDR; standard for biological
  ##   data where missing a real outlier costs more than flagging a
  ##   marginal one.
  ## - "gap": per-channel data-driven break via the package's broken-
  ##   stick + tail-decay inflection helper.  Useful for diagnostic
  ##   exploration but mis-calibrates on the half-normal / Rayleigh
  ##   null (the broken-stick model finds breaks in the natural decay
  ##   of thin-tailed nulls, over-flagging clean tracks).  Kept as
  ##   opt-in.
  ##
  ## Under variance = "dbbmm" only the chisq channel exists (no per-
  ## axis decomposition); per-axis thresholds stay NA.
  z_thr_para  <- NA_real_
  z_thr_orth  <- NA_real_
  z_thr_chisq <- NA_real_

  if (z_threshold_method == "gap") {
    if (!is.null(z_threshold)) {
      rlang::abort(paste0(
        "`z_threshold` numeric override is incompatible with ",
        "z_threshold_method = 'gap'; the gap method is data-driven. ",
        "Either drop z_threshold or switch to a parametric method."),
        class = "move2utils_mt_flag_outliers_dbgb_z_threshold_with_gap")
    }
    ## NULL `z_gap_threshold` defers to .z_gap_threshold's leaf formal
    ## (single source of truth).  Resolve once for messaging then
    ## forward the resolved value to the three leaf calls.
    z_gt_eff <- if (is.null(z_gap_threshold))
                  formals(.z_gap_threshold)$threshold
                else z_gap_threshold
    if (variance == "dbgb") {
      z_thr_para  <- .z_gap_threshold(bridge$Z_para,  z_gt_eff)
      z_thr_orth  <- .z_gap_threshold(bridge$Z_orth,  z_gt_eff)
    }
    z_thr_chisq <- .z_gap_threshold(bridge$Z_chisq, z_gt_eff)
    say(sprintf(
      "Z thresholds (gap, threshold=%g): para=%s orth=%s chisq=%s",
      z_gt_eff,
      if (is.finite(z_thr_para))  sprintf("%.3f", z_thr_para)  else "Inf",
      if (is.finite(z_thr_orth))  sprintf("%.3f", z_thr_orth)  else "Inf",
      if (is.finite(z_thr_chisq)) sprintf("%.3f", z_thr_chisq) else "Inf"))
  } else {
    fwer_total <- 0.05
    fdr_total  <- 0.05
    ## Numeric override applies to the chisq channel only.  Per-axis
    ## channels stay parametric (or NA under variance = "dbbmm").
    if (!is.null(z_threshold)) {
      if (!is.numeric(z_threshold) || length(z_threshold) != 1L ||
          is.na(z_threshold) || z_threshold <= 0)
        rlang::abort("`z_threshold` must be a positive scalar or NULL.",
                     class = "move2utils_mt_flag_outliers_dbgb_bad_z_threshold")
      z_thr_chisq <- z_threshold
    } else {
      if (z_threshold_method == "bh_fdr") {
        z_thr_chisq <- .bh_fdr_z(bridge$Z_chisq, fdr = fdr_total / 2)
      } else {
        z_thr_chisq <- .bonferroni_z(ns, fwer = fwer_total / 2)
      }
    }
    if (variance == "dbgb") {
      if (z_threshold_method == "bh_fdr") {
        z_thr_para <- .bh_fdr_z_normal(bridge$Z_para, fdr = fdr_total / 4)
        z_thr_orth <- .bh_fdr_z_normal(bridge$Z_orth, fdr = fdr_total / 4)
      } else {
        z_thr_para <- .bonferroni_z_normal(ns, fwer = fwer_total / 4)
        z_thr_orth <- .bonferroni_z_normal(ns, fwer = fwer_total / 4)
      }
    }
    say(sprintf(
      "Z thresholds (%s envelope, joint FWER<=%g): para=%s orth=%s chisq=%s",
      z_threshold_method, fwer_total,
      if (is.finite(z_thr_para))  sprintf("%.3f", z_thr_para)  else "Inf/NA",
      if (is.finite(z_thr_orth))  sprintf("%.3f", z_thr_orth)  else "Inf/NA",
      if (is.finite(z_thr_chisq)) sprintf("%.3f", z_thr_chisq) else "Inf"))
  }

  ## Per-channel exceedance.
  fl_para <- is.finite(bridge$Z_para)  & bridge$Z_para  > z_thr_para
  fl_orth <- is.finite(bridge$Z_orth)  & bridge$Z_orth  > z_thr_orth
  fl_chisq <- is.finite(bridge$Z_chisq) & bridge$Z_chisq > z_thr_chisq
  ## NA-vs-FALSE hygiene: treat NA channels as not-flagged.
  fl_para[is.na(fl_para)] <- FALSE
  fl_orth[is.na(fl_orth)] <- FALSE
  fl_chisq[is.na(fl_chisq)] <- FALSE
  flagged_by_z <- fl_para | fl_orth | fl_chisq

  ## Diagnostic label per fix: which channel(s) fired plus dominant
  ## axis when only chisq fires.  This refinement matters under
  ## parametric calibration where per-axis tests rarely fire (the
  ## per-axis Bonferroni threshold sits well above the natural Z bulk
  ## edge): the chisq channel is doing most of the work, so a fix
  ## flagged purely by chisq still carries directional information
  ## via the (Z_para, Z_orth) ratio.  A chisq-only fix with
  ## |Z_para| > 2|Z_orth| reads as "isotropic_para_dom" -- the
  ## anomaly was driven by the along-axis component even though the
  ## per-axis test wasn't extreme enough to fire on its own.
  abs_zp <- abs(bridge$Z_para)
  abs_zo <- abs(bridge$Z_orth)
  abs_zp[is.na(abs_zp)] <- 0
  abs_zo[is.na(abs_zo)] <- 0
  para_dom <- abs_zp > 2 * abs_zo
  orth_dom <- abs_zo > 2 * abs_zp
  bridge_class_surv <- rep("none", ns)
  bridge_class_surv[ fl_para & !fl_orth & !fl_chisq] <- "para"
  bridge_class_surv[!fl_para &  fl_orth & !fl_chisq] <- "orth"
  ## chisq-only: refine by axis dominance ratio
  chisq_only <- !fl_para & !fl_orth & fl_chisq
  bridge_class_surv[chisq_only & para_dom] <- "isotropic_para_dom"
  bridge_class_surv[chisq_only & orth_dom] <- "isotropic_orth_dom"
  bridge_class_surv[chisq_only & !para_dom & !orth_dom] <- "isotropic"
  bridge_class_surv[ fl_para &  fl_orth & !fl_chisq] <- "para+orth"
  bridge_class_surv[ fl_para & !fl_orth &  fl_chisq] <- "para+chisq"
  bridge_class_surv[!fl_para &  fl_orth &  fl_chisq] <- "orth+chisq"
  bridge_class_surv[ fl_para &  fl_orth &  fl_chisq] <- "all"

  ## ---- resolve residual_max (absolute-residual axis) ----
  ## NULL  = auto via gap detection on log(|residual|).  Catches
  ##         geometrically extreme residuals that dBBMM-Z would absorb
  ##         into local sigma^2.
  ## Inf   = disabled (Z-only criterion).
  ## numeric = hard cap.
  if (is.null(residual_max)) {
    residual_max <- .auto_residual_max(
      bridge$residual,
      entropy_threshold = residual_entropy_threshold,
      gap_threshold     = residual_gap_threshold,
      dip_alpha         = residual_dip_alpha)
    say(sprintf("residual_max (auto, gap detection on log(residual)): %s",
                if (is.finite(residual_max))
                  sprintf("%.1f m", residual_max) else "no break"))
  }
  if (!is.numeric(residual_max) || length(residual_max) != 1L) {
    rlang::abort("`residual_max` must be a positive scalar, Inf, or NULL.",
                 class = "move2utils_mt_flag_outliers_dbgb_bad_residual_max")
  }
  flagged_by_residual <- is.finite(bridge$residual) &
                          is.finite(residual_max) &
                          bridge$residual > residual_max
  flagged_in_surv <- flagged_by_z | flagged_by_residual

  ## Fold residual-only flags into the diagnostic taxonomy: a fix
  ## flagged exclusively by the residual gate (Z channels all below
  ## their thresholds) gets class "residual".  This is the geometric-
  ## impossibility cap firing on extreme residuals that dBBMM-Z would
  ## absorb into local sigma^2.
  bridge_class_surv[flagged_by_residual & !flagged_by_z] <- "residual"

  ## ---- assemble full-length result columns ----
  bridge_residual    <- rep(NA_real_, n_total)
  bridge_width       <- rep(NA_real_, n_total)
  sigma2_motion_para <- rep(NA_real_, n_total)
  sigma2_motion_orth <- rep(NA_real_, n_total)
  sigma2_motion      <- rep(NA_real_, n_total)
  bridge_z           <- rep(NA_real_, n_total)
  bridge_z_para      <- rep(NA_real_, n_total)
  bridge_z_orth      <- rep(NA_real_, n_total)
  bridge_z_chisq     <- rep(NA_real_, n_total)
  bridge_z_class     <- rep(NA_character_, n_total)
  bridge_iter        <- rep(NA_integer_, n_total)
  is_outlier         <- is_peeled

  bridge_residual[surv_idx] <- bridge$residual
  bridge_width[surv_idx]    <- sqrt(bridge$width2)
  if (variance == "dbgb") {
    sigma2_motion_para[surv_idx] <- mv$para
    sigma2_motion_orth[surv_idx] <- mv$orth
    bridge_z_para[surv_idx]      <- bridge$Z_para
    bridge_z_orth[surv_idx]      <- bridge$Z_orth
    bridge_z_chisq[surv_idx]     <- bridge$Z_chisq
  } else {
    sigma2_motion[surv_idx]   <- mv
    bridge_z_chisq[surv_idx]  <- bridge$Z_chisq
  }
  bridge_z[surv_idx]            <- bridge$Z
  bridge_z_class[surv_idx]      <- bridge_class_surv
  is_outlier[surv_idx[flagged_in_surv]] <- TRUE
  bridge_iter[is_peeled] <- 0L
  bridge_iter[surv_idx[flagged_in_surv]] <- 1L

  x$bridge_residual  <- bridge_residual
  x$bridge_width     <- bridge_width
  if (variance == "dbgb") {
    x$sigma2_motion_para <- sigma2_motion_para
    x$sigma2_motion_orth <- sigma2_motion_orth
    x$bridge_z_para      <- bridge_z_para
    x$bridge_z_orth      <- bridge_z_orth
  } else {
    x$sigma2_motion    <- sigma2_motion
  }
  x$bridge_z         <- bridge_z
  x$bridge_z_chisq   <- bridge_z_chisq
  x$bridge_z_class   <- bridge_z_class
  x$bridge_iteration <- bridge_iter
  x$is_outlier       <- is_outlier

  ## Per-channel thresholds: para / orth are NA in non-envelope modes
  ## and under variance = "dbbmm".  z_threshold (singular) preserved
  ## as the binding chisq threshold for back-compat with prior 0.2.0
  ## attribute consumers.
  attr(x, "z_thresholds") <- list(para  = z_thr_para,
                                   orth  = z_thr_orth,
                                   chisq = z_thr_chisq)
  attr(x, "z_threshold")  <- z_thr_chisq
  attr(x, "z_threshold_method") <- z_threshold_method
  attr(x, "residual_max") <- residual_max
  attr(x, "window_size")  <- window_size
  attr(x, "location_error")      <- le_full

  say(sprintf("=== mt_flag_outliers_dbgb: %d flagged (%.3f%% of %d); peel=%d, by-channel para=%d orth=%d chisq=%d, residual=%d ===",
              sum(is_outlier), 100 * sum(is_outlier) / n_total, n_total,
              sum(is_peeled),
              sum(fl_para), sum(fl_orth), sum(fl_chisq),
              sum(flagged_by_residual)))

  if (orig_was_longlat) {
    cols <- if (variance == "dbgb")
      c("bridge_residual","bridge_width",
        "sigma2_motion_para","sigma2_motion_orth",
        "bridge_z","bridge_z_chisq","bridge_z_class",
        "bridge_z_para","bridge_z_orth",
        "bridge_iteration","is_outlier")
    else
      c("bridge_residual","bridge_width","sigma2_motion",
        "bridge_z","bridge_z_chisq","bridge_z_class",
        "bridge_iteration","is_outlier")
    for (col in cols) orig_x[[col]] <- x[[col]]
    attr(orig_x, "z_thresholds")       <- attr(x, "z_thresholds")
    attr(orig_x, "z_threshold")        <- attr(x, "z_threshold")
    attr(orig_x, "z_threshold_method") <- attr(x, "z_threshold_method")
    attr(orig_x, "residual_max")       <- attr(x, "residual_max")
    attr(orig_x, "window_size")        <- attr(x, "window_size")
    attr(orig_x, "location_error")     <- attr(x, "location_error")
    x <- orig_x
  }

  if (plot) .plot_bridge_dbbmm(x, z_thr_chisq)
  if (remove) x <- x[!x$is_outlier, ]
  x
}


# ---- helpers ---------------------------------------------------------------

## Bonferroni Z threshold for FWER alpha across n active fixes.
## Under the bridge null, Z ~ Rayleigh(1) so P(Z > z) = exp(-z^2/2).
## Per-fix alpha = fwer / n -> z = sqrt(-2 log(fwer / n)).
## @keywords internal
.bonferroni_z <- function(n, fwer = 0.05) {
  n <- max(1L, as.integer(n))
  sqrt(-2 * log(fwer / n))
}

## Benjamini-Hochberg FDR Z threshold given a vector of Z scores.
##
## Under the bridge null Z ~ Rayleigh(1), so the per-fix p-value is
## p_i = exp(-Z_i^2 / 2).  BH at level alpha rejects all p-values
## below the line k/n * alpha (sorted ascending).  The Z threshold
## is the smallest Z whose p-value gets rejected.
##
## FDR controls the EXPECTED FRACTION of false positives among the
## flagged set, which is the right operating point when missing a
## real outlier is more costly than flagging a marginal one --
## standard for biological data where FN cost dominates.  More
## sensitive than Bonferroni FWER at the same nominal alpha.
##
## Returns Inf when no Z exceeds any BH cut (no flags justified).
##
## @keywords internal
.bh_fdr_z <- function(Z, fdr = 0.05) {
  Z <- Z[is.finite(Z) & Z > 0]
  n <- length(Z)
  if (n < 2L) return(Inf)
  ## p-values under Rayleigh(1)
  p <- exp(-Z^2 / 2)
  p_sorted <- sort(p)
  bh_line <- (seq_len(n) / n) * fdr
  rejected <- p_sorted <= bh_line
  if (!any(rejected)) return(Inf)
  k <- max(which(rejected))
  ## Z threshold = the Z at the k-th smallest p
  sqrt(-2 * log(p_sorted[k]))
}

## Data-driven gap threshold on a Z-channel using the package's
## broken-stick + tail-decay inflection detector.
##
## The gap detector (.gap_threshold_lower) finds the lower-tail break
## in a numeric vector where outliers are at the LOW end.  Z scores
## have outliers at the HIGH end, so we feed -log(Z) -- mirrors the
## convention used by mt_flag_outliers_bridge on bridge_eta.  The log
## transform compresses the heavy upper tail of contaminated Z so the
## broken-stick / inflection logic operates on a well-behaved signal.
##
## Returns the Z threshold (in raw Z units) above which fixes are
## flagged on that channel; Inf when the detector finds no defensible
## break (clean data).
##
## @keywords internal
.z_gap_threshold <- function(Z, threshold = 3) {
  Z <- Z[is.finite(Z) & Z > 0]
  if (length(Z) < 10L) return(Inf)
  log_scores <- -log(Z)
  gap <- .gap_threshold_lower(log_scores, threshold = threshold)
  if (is.na(gap$break_value)) return(Inf)
  ## break_value is in -log(Z) space; convert back to Z.
  exp(-gap$break_value)
}


## Bonferroni axis-Z threshold for FWER alpha across n axis tests.
## Per-axis Z is |Normal(0,1)| under the dBGB null (each component of
## the residual normalised by its own axis sigma).  Two-sided per-axis
## p = 2 * (1 - Phi(|Z|)) gives the half-normal tail; Bonferroni cut
## per-axis alpha = fwer / n -> z = qnorm(1 - fwer / (2n)).
## n here is the total axis test count (= 2 * fixes for the union mode).
## @keywords internal
.bonferroni_z_normal <- function(n, fwer = 0.05) {
  n <- max(1L, as.integer(n))
  stats::qnorm(1 - fwer / (2 * n))
}

## Benjamini-Hochberg FDR axis-Z threshold given pooled axis Z-scores.
##
## Z is the concatenation c(Z_para, Z_orth).  Under the dBGB null each
## entry is |Normal(0,1)|; the per-axis two-sided p-value is
## p_i = 2 * (1 - Phi(|Z_i|)).  BH at level alpha rejects p_i below
## (k/n) * alpha (sorted ascending); the Z threshold is the |Normal|
## quantile at the k-th p.  Returns Inf when no axis test exceeds any
## BH cut.
##
## @keywords internal
.bh_fdr_z_normal <- function(Z, fdr = 0.05) {
  Z <- Z[is.finite(Z) & Z > 0]
  n <- length(Z)
  if (n < 2L) return(Inf)
  ## two-sided per-axis p under |N(0,1)|
  p <- 2 * stats::pnorm(-abs(Z))
  p_sorted <- sort(p)
  bh_line <- (seq_len(n) / n) * fdr
  rejected <- p_sorted <= bh_line
  if (!any(rejected)) return(Inf)
  k <- max(which(rejected))
  stats::qnorm(1 - p_sorted[k] / 2)
}

## Auto window_size from ACF of v^2 (motion intensity).
##
## v^2 = step^2/dt is the per-step motion-intensity score.  Its
## autocorrelation function decays from 1 at lag 0 toward 0 at lags
## beyond the typical motion-coherence timescale.  Because v^2 is
## heavy-tailed, lag-1 ACF is already substantially below 1 even on
## smooth tracks (the heavy tail decorrelates immediately at lag 1
## while the bulk persists).  An ACF threshold of 0.5 is therefore
## too generous; 0.2 captures the timescale at which the bulk has
## decorrelated.
##
## window = c_factor * lag_decorrelation, capped at min(max_window,
## n/8) to keep dBBMM tractable.  The c_factor of 4 captures multiple
## decorrelated motion episodes per window so BIC has data to find
## breakpoints.
##
## @keywords internal
.auto_window_size_acf <- function(x,
                                    c_factor       = 4,
                                    acf_threshold  = 0.2,
                                    min_window     = 21L,
                                    max_window     = 201L,
                                    default_window = 49L) {
  cc <- sf::st_coordinates(x)
  ts <- as.numeric(move2::mt_time(x), units = "secs")
  n  <- nrow(cc)
  if (n < 30L) return(default_window)
  step_m <- sqrt(diff(cc[, 1L])^2 + diff(cc[, 2L])^2)
  dt_s   <- diff(ts)
  v2 <- step_m^2 / dt_s
  v2 <- v2[is.finite(v2) & v2 > 0]
  if (length(v2) < 30L) return(default_window)
  max_lag <- min(200L, length(v2) %/% 4L)
  acf_out <- stats::acf(v2, lag.max = max_lag, plot = FALSE,
                          na.action = stats::na.pass)
  acf_vals <- as.numeric(acf_out$acf)[-1L]    # lag 0 is always 1
  lag_decorr <- which(acf_vals < acf_threshold)[1L]
  if (is.na(lag_decorr)) lag_decorr <- max_lag
  win <- as.integer(c_factor * lag_decorr)
  win <- max(min_window, min(win, max_window))
  if (win %% 2L == 0L) win <- win - 1L
  win
}


## Auto residual_max via gap detection on log(|residual|).
##
## Returns the threshold residual magnitude in metres above which
## fixes are flagged regardless of state-aware Z.  Uses the package's
## existing entropy-valley + dip-validated broken-stick detector on
## -log(residual) (LOWER tail in -log space = UPPER tail in residual
## space).
##
## When the residual distribution is multi-scale (outliers themselves
## span orders of magnitude), the detector can place the break inside
## the outlier set rather than at the bulk-tail boundary, missing the
## smaller-but-still-real outlier subset.  The user can tighten the
## gap_threshold parameter to push the break toward the bulk for
## tracks where missing any real outlier is unacceptable.
##
## @keywords internal
.auto_residual_max <- function(residual,
                                 entropy_threshold = NULL,
                                 gap_threshold     = NULL,
                                 dip_alpha         = 0.05) {
  ## NULL defers to leaf formals; non-NULL overrides.  See
  ## audits/2026-05-25-parameter-propagation/findings.md §1.2.
  r <- residual[is.finite(residual) & residual > 0]
  if (length(r) < 30L) return(Inf)
  br <- .entropy_or_dip_gap_threshold_lower(
    -log(r),
    entropy_threshold = entropy_threshold,
    gap_threshold     = gap_threshold,
    dip_alpha         = dip_alpha)
  if (is.na(br$break_value)) return(Inf)
  exp(-br$break_value)
}


## Compute bridge residual, width and Z given sigma^2 series.
## @keywords internal
.compute_bridge_z <- function(x, sigma2_per_min, location_error) {
  cc <- sf::st_coordinates(x)
  ts <- as.numeric(move2::mt_time(x), units = "secs")
  n  <- nrow(cc)
  dt1 <- c(NA_real_, diff(ts))
  dt2 <- c(diff(ts), NA_real_)
  prev_ <- c(NA_integer_, 1:(n - 1L))
  next_ <- c(2:n, NA_integer_)
  mean_x <- (cc[prev_, 1L] * dt2 + cc[next_, 1L] * dt1) / (dt1 + dt2)
  mean_y <- (cc[prev_, 2L] * dt2 + cc[next_, 2L] * dt1) / (dt1 + dt2)
  res <- sqrt((cc[, 1L] - mean_x)^2 + (cc[, 2L] - mean_y)^2)
  width2 <- (dt1 * dt2) / (dt1 + dt2)
  sigma2_per_s <- sigma2_per_min / 60
  expected_var <- sigma2_per_s * width2 + 2 * location_error^2
  Z <- res / sqrt(expected_var)
  list(residual = res, width2 = width2, Z = Z)
}


## Directional analogue of .compute_bridge_z.
##
## sigma2_per_min is a 2-column data.frame (para, orth) per fix from
## mt_motion_variance.mt_dbgb_variance(); each column is in m^2/min,
## the dBBMM/dBGB convention.  Decomposes the residual onto the local
## travel axis (mu -> next_fix) using the same convention as
## .compute_bridge_residuals_dbgb in mt_flag_outliers_bridge.R and
## .delta_para_orth in dbgb_variance.R, then normalises each component
## by its own per-axis sigma:
##
##   Z_para = |r_para| / sqrt(sigma2_para_per_s * width^2 + 2*loc_err^2)
##   Z_orth = |r_orth| / sqrt(sigma2_orth_per_s * width^2 + 2*loc_err^2)
##
## Each Z_axis ~ |Normal(0,1)| under the bridge null with axis-sigma.
## sqrt(Z_para^2 + Z_orth^2) ~ Rayleigh(1) (chi^2_2 magnitude), the
## same null distribution as the dbbmm scalar Z.  Returns also the
## isotropic Euclidean residual so the residual_max gate (geometric-
## impossibility cap) stays scalar.
##
## When the bridge axis is degenerate (mu == next_fix), the residual
## is split equally (sqrt(2) factor) between para and orth -- mirrors
## .compute_bridge_residuals_dbgb's "no axis direction" fallback.
## @keywords internal
.compute_bridge_z_dbgb <- function(x, sigma2_per_min, location_error) {
  if (!is.data.frame(sigma2_per_min) ||
      !all(c("para", "orth") %in% names(sigma2_per_min)))
    rlang::abort(
      "sigma2_per_min must be a data.frame with columns 'para' and 'orth'.",
      class = "move2utils_mt_flag_outliers_dbgb_bad_sigma2_per_min")
  cc <- sf::st_coordinates(x)
  ts <- as.numeric(move2::mt_time(x), units = "secs")
  n  <- nrow(cc)
  dt1 <- c(NA_real_, diff(ts))
  dt2 <- c(diff(ts), NA_real_)
  prev_ <- c(NA_integer_, 1:(n - 1L))
  next_ <- c(2:n, NA_integer_)
  mean_x <- (cc[prev_, 1L] * dt2 + cc[next_, 1L] * dt1) / (dt1 + dt2)
  mean_y <- (cc[prev_, 2L] * dt2 + cc[next_, 2L] * dt1) / (dt1 + dt2)

  rx <- cc[, 1L] - mean_x
  ry <- cc[, 2L] - mean_y
  rn <- sqrt(rx * rx + ry * ry)

  ## decomposition axis: unit vector mu -> next fix.  Same convention
  ## as .delta_para_orth in dbgb_variance.R and the dBGB residual
  ## helper in mt_flag_outliers_bridge.R, so per-axis sigmas align
  ## with the variance estimator.
  ax_x <- cc[next_, 1L] - mean_x
  ax_y <- cc[next_, 2L] - mean_y
  axn  <- sqrt(ax_x * ax_x + ax_y * ax_y)

  has_axis <- is.finite(axn) & axn > 0
  rp <- rep(NA_real_, n)
  ro <- rep(NA_real_, n)
  rp[has_axis] <- abs(rx[has_axis] * ax_x[has_axis] +
                       ry[has_axis] * ax_y[has_axis]) / axn[has_axis]
  ro[has_axis] <- sqrt(pmax(rn[has_axis] ^ 2 - rp[has_axis] ^ 2, 0))
  ## degenerate axis: split equally between para and orth.
  no_axis <- is.finite(rn) & !has_axis
  rp[no_axis] <- rn[no_axis] / sqrt(2)
  ro[no_axis] <- rn[no_axis] / sqrt(2)

  width2 <- (dt1 * dt2) / (dt1 + dt2)
  sigma2_para_per_s <- sigma2_per_min$para / 60
  sigma2_orth_per_s <- sigma2_per_min$orth / 60
  expected_var_para <- sigma2_para_per_s * width2 + 2 * location_error^2
  expected_var_orth <- sigma2_orth_per_s * width2 + 2 * location_error^2
  Z_para <- rp / sqrt(expected_var_para)
  Z_orth <- ro / sqrt(expected_var_orth)
  ## Z left NA at the call site; the caller fills it from
  ## (Z_para, Z_orth) according to z_aggregate.
  list(residual = rn, width2 = width2,
       Z = rep(NA_real_, n),
       Z_para = Z_para, Z_orth = Z_orth)
}


## Diagnostic plot for the dBBMM-Z primitive.
## @keywords internal
.plot_bridge_dbbmm <- function(x, z_threshold) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  Z <- x$bridge_z
  Z_pos <- Z[is.finite(Z) & Z > 0]
  if (length(Z_pos) >= 10L) {
    d <- stats::density(log(Z_pos), n = 512)
    graphics::plot(d$x, d$y, type = "l", col = "grey20", lwd = 1.4,
                   xlab = "log(Z)", ylab = "density",
                   main = sprintf("dBBMM-Z distribution (n=%d)",
                                  length(Z_pos)))
    graphics::abline(v = log(z_threshold), col = "firebrick",
                     lwd = 2, lty = 2)
    graphics::legend("topright",
                     legend = sprintf("Z_thr = %.2f", z_threshold),
                     col = "firebrick", lty = 2, lwd = 2, bty = "n",
                     cex = 0.85)
  } else {
    graphics::plot.new()
    graphics::title("dBBMM-Z distribution")
  }

  cc <- sf::st_coordinates(x)
  graphics::plot(cc, type = "l", col = "grey70", lwd = 0.3, asp = 1,
                 xlab = "", ylab = "",
                 main = sprintf("flagged: %d (%.3f%%)",
                                sum(x$is_outlier, na.rm = TRUE),
                                100 * mean(x$is_outlier, na.rm = TRUE)))
  flagged <- which(x$is_outlier)
  if (length(flagged)) {
    graphics::points(cc[flagged, , drop = FALSE],
                     col = grDevices::adjustcolor("red", 0.6),
                     pch = 20, cex = 0.7)
  }
  invisible(NULL)
}

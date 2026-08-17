#' Flag or remove outliers in movement data based on joint movement probabilities
#'
#' Detects outliers in a \code{move2} object by computing joint probabilities
#' from the empirical distributions of step lengths, turning angles, and their
#' consecutive changes. Locations that fall in low-probability regions of this
#' joint space are flagged as potential outliers.
#'
#' @details
#' The method works by building three probability components for each location.
#' The three components are: (1) step-turn probability, from a 2D histogram
#' of step lengths and turning angles with circular wrapping and bilinear
#' interpolation; (2) delta-step probability, from the kernel density of
#' changes in step length between consecutive steps; and (3) delta-turn
#' probability, from the kernel density of changes in turning angle between
#' consecutive steps.
#'
#' The joint probability is the product of all three. Locations whose joint
#' probability falls below the specified percentile threshold are flagged as
#' outliers.
#'
#' When the input contains multiple individuals, each is processed separately
#' by default. To build a pooled reference distribution from all individuals,
#' pass \code{reference = x}. To use an external clean track as the reference,
#' pass it via \code{reference}.
#'
#' Two methods are available for computing the step-turn probability:
#' \code{"histogram"} (default) uses a 2D histogram with circular wrapping
#' and bilinear interpolation; \code{"copula"} fits parametric marginal
#' distributions (Weibull for step lengths, von Mises for turning angles)
#' and uses the product of marginal densities. The copula method is faster
#' and can work better with small samples.
#'
#' When \code{iterations > 1}, the detection runs iteratively: after each
#' pass, flagged outliers are masked (removed from the track) and movement
#' metrics are recomputed on the cleaned track. This allows detection of
#' consecutive outliers, because removing the first outlier reveals the
#' true step from the last good location to the next good location.
#' Iteration stops after the specified number of passes or when no new
#' outliers are found, whichever comes first. All flags are mapped back
#' to the original object.
#'
#' When \code{quality_columns} is provided, each quality function maps a
#' raw data column to a \[0,1\] quality score. The product of all quality
#' scores multiplies the movement probability before thresholding:
#' \deqn{P_{final} = P_{movement} \times \prod quality\_weights}
#'
#' When \code{time_normalize = TRUE}, the method uses speed
#' (step_length / time_lag) and angular velocity (turning_angle / time_lag)
#' instead of raw step lengths and turning angles. This makes the method
#' time-aware, which is recommended for irregularly sampled data.
#'
#' @param x A \code{move2} object. Must contain at least 3 non-empty
#'   locations.  Either lon/lat or projected; projected input is
#'   transformed internally to WGS84 lon/lat for the turning-angle
#'   computation (see \code{move2::mt_azimuth}) and the result is
#'   returned in the original CRS.
#' @param threshold Numeric, or \code{NULL} for automatic defaults.
#'   Interpretation depends on \code{threshold_type}:
#'   With \code{"gap"} (default): how many times the local spacing a gap in
#'   the sorted log-probabilities must exceed to be considered a natural break.
#'   Default is 3. Higher values are more conservative.
#'   With \code{"entropy"}: the maximum allowed valley-to-peak density ratio
#'   in the KDE of log-probabilities; a valley deeper than this declares
#'   an outlier regime. Default is 0.3 (unified package-wide entropy
#'   default; see \code{\link{mt_flag_outliers_bridge}} for the same
#'   value's rationale). Higher values admit shallower valleys (more
#'   sensitive); lower values demand a deeper split.
#'   With \code{"significance"}: the significance level for flagging based on
#'   robust z-scores of log-probabilities. Default is 0.001.
#'   With \code{"percentile"}: the bottom fraction to flag (e.g. 0.001 flags
#'   the bottom 0.1 percent). Default is 0.001.
#'   When \code{NULL}, the default for the chosen \code{threshold_type} is used.
#' @param threshold_type Character. One of \code{"gap"} (default),
#'   \code{"entropy"}, \code{"significance"}, or \code{"percentile"}.
#'   The \strong{gap} method uses a two-stage approach: (1) a broken stick
#'   null model (MacArthur 1957) screens for candidate breaks by comparing
#'   observed gap sizes in the sorted log-probabilities to their expected
#'   sizes under a single continuous distribution; (2) a tail-decay
#'   inflection analysis tracks how the tail shortens as points are removed
#'   from the extreme left — genuine outliers cause steep drops, while
#'   entering the bulk distribution produces small, linear changes. The
#'   inflection point where the second derivative of the tail-length curve
#'   drops to the noise floor determines the natural boundary. No
#'   distributional assumptions are made, and clean data produces no or
#'   very few outliers. Default threshold is 3.
#'   The \strong{entropy} method estimates the density of log-probabilities
#'   via KDE and searches for the deepest local minimum (valley) below the
#'   main mode; a fix is flagged when its log-probability sits below that
#'   valley. The default \code{threshold = 0.3} is the unified package-
#'   wide maximum valley-to-peak density ratio (see
#'   \code{\link{mt_flag_outliers_bridge}} for the same value's
#'   rationale and the Raven-sweep validation). Returns no outliers
#'   on clean unimodal data.
#'   The \strong{significance} method uses robust z-scores (median + MAD) on
#'   log-probabilities, assuming approximate normality. It can over-flag in
#'   left-skewed distributions.
#'   The \strong{percentile} method always flags the bottom fraction
#'   regardless of data quality and does not converge during iteration.
#' @param prob_type Character. Which probability to use for outlier detection.
#'   One of \code{"joint"} (default), \code{"step_turn"}, \code{"delta_step"},
#'   \code{"delta_turn"}, or \code{"custom"} (product of step_turn and
#'   delta_step only).
#' @param remove Logical. If \code{TRUE}, return the object with outliers
#'   removed. If \code{FALSE} (default), return the original object with
#'   outlier flags and probabilities added as columns.
#' @param plot Logical. If \code{TRUE} (default), create a two-panel diagnostic
#'   plot showing probability-coloured locations and flagged outliers.
#' @param autodiff_alpha Exponent on the auto-difference terms in the
#'   joint probability (Equation 1 in the paper). Controls how much
#'   weight the persistence components (delta-speed, delta-angular
#'   velocity) carry relative to the step--turn component.
#'
#'   The default \code{"acf"} derives alpha from the data: the lag-1
#'   autocorrelation of speed (\eqn{r_v}) and angular velocity
#'   (\eqn{r_\omega}) are computed, and \eqn{\alpha = \sqrt{r_v \cdot
#'   r_\omega}}.  High autocorrelation means persistence is informative
#'   and auto-differences should carry weight; low autocorrelation means
#'   they are noise and should be downweighted.
#'
#'   Numeric values override the ACF estimate: \code{0} ignores
#'   auto-differences, \code{0.5} applies a geometric mean, \code{1.0}
#'   gives the simple product.  The string \code{"auto"} selects alpha
#'   by maximising the trimmed mean of log-probabilities (legacy mode).
#' @param method Character. Method for computing step-turn probabilities.
#'   One of \code{"histogram"} (default, 2D histogram with bilinear
#'   interpolation) or \code{"copula"} (parametric marginal distributions:
#'   Weibull for step lengths, von Mises for turning angles). The copula
#'   method requires the \code{circular} and \code{MASS} packages.
#' @param step_transform Character. Controls the axis on which the 2D
#'   turn/step histogram is built when \code{method = "histogram"}.
#'   \code{"none"} (default) uses the raw step length and preserves the
#'   joint turn/step structure that the synthetic benchmark was
#'   validated on. \code{"log"} applies \code{log(1 + step)} to the
#'   y-axis; useful only for tracks with a pathologically heavy
#'   step-length tail (e.g. teleport-class GPS errors spanning several
#'   orders of magnitude) where raw-scale binning collapses real
#'   movement into a single row of the histogram. The log transform is
#'   a monotonic, invertible change of variable, not a distributional
#'   assumption — the density is still estimated non-parametrically
#'   from the transformed data. Note: log-transform can hide
#'   physiologically-plausible joint outliers (wrong-turn / impossible
#'   acceleration at reasonable step lengths), so it is not the
#'   default. If you pass a track whose
#'   \code{diff(range(step))/IQR(step)} is very large, the function
#'   emits a suggestion to try \code{"log"}. Ignored when
#'   \code{method = "copula"}, where the parametric Weibull marginal
#'   handles the tail on its own.
#' @param step_floor Numeric, non-negative. Minimum absolute step length
#'   (metres) for a rate-flagged fix to actually be flagged as an
#'   outlier. Default \code{0} (disabled); the pre-2026 rate-only
#'   behaviour. Set to a positive value (commonly 5--25 m, or your
#'   device's nominal accuracy) to add a two-axis criterion: a fix
#'   is flagged only where both the joint-probability threshold is
#'   tripped AND the absolute step length exceeds this floor.
#'   Recommended for real-world data with burst-mode sampling, where
#'   a few-metre GPS jitter over a 1-second dt produces a "high
#'   speed" that is not a physical outlier. Not on by default so
#'   existing synthetic-ground-truth benchmarks (which can include
#'   sub-noise displacement outliers by construction) continue to
#'   pass. Applied after thresholding, before iterative refinement.
#' @param reference Optional \code{move2} object from which to build the
#'   probability surfaces. If \code{NULL} (default), each individual's own
#'   data are used.  Supply a longer or cleaner track to improve
#'   distributions for short or contaminated tracks.  To pool all
#'   individuals into a single reference distribution, pass
#'   \code{reference = x}.  Mutually exclusive with \code{pool_by}.
#' @param pool_by Optional character vector of length 1 or 2 naming
#'   column(s) in \code{mt_track_data(x)}.  Length 1: single column
#'   used as the pool source (e.g. \code{"individual_id"}).
#'   Length 2: \code{c(outer, inner)} where \code{outer} names the
#'   fit-source column -- the union of its events supplies one
#'   \code{(step, turn, delta_step, delta_turn, gaps)} reference
#'   distribution per outer group, injected into each member track's
#'   per-track dispatch.  Length 2 requires strict nesting: every
#'   distinct \code{inner} value must map to exactly one
#'   \code{outer} value.  Length \eqn{> 2} is rejected.  Note: for
#'   the probability primitive the pool is \emph{integrated} (acts
#'   through the per-track dispatcher) rather than post-hoc, so
#'   \code{inner} has no role here -- it is validated for
#'   consistency with the orchestrator but does not affect the
#'   prob primitive's flagging beyond what \code{outer} drives.
#'   \code{NULL} (default) preserves per-track behaviour.  NA
#'   values in the named column(s) cause those tracks to fall back
#'   to per-track processing with a warning.  Mutually exclusive
#'   with \code{reference}.
#' @param drop_na Logical. If \code{TRUE}, also remove locations where the
#'   probability could not be calculated (e.g. first/last locations). Default
#'   is \code{FALSE}.
#' @param iterations Integer. Number of iterative refinement passes. Default
#'   is \code{1} (no iteration). When greater than 1, after each pass flagged
#'   outliers are masked and movement metrics are recomputed on the cleaned
#'   track. This helps detect consecutive outliers. Iteration stops early if
#'   no new outliers are found.
#' @param quality_columns Named list of functions, or \code{NULL} (default).
#'   Each name must be a column in \code{x}, and each function maps the raw
#'   column values to a \[0,1\] quality score. The product of all quality
#'   scores multiplies the joint probability before thresholding.
#'   Example:
#'   \preformatted{quality_columns = list(
#'     "gps.satellite.count" = function(s) pnorm(s, mean=7, sd=2),
#'     "gps.hdop" = function(h) 1 - pnorm(h, mean=3, sd=1.5)
#'   )}
#' @param time_normalize Logical. If \code{TRUE} (default), use speed
#'   (step_length / time_lag) and angular velocity (turning_angle / time_lag)
#'   instead of raw step lengths and turning angles. This makes the method
#'   time-aware: the same displacement over different time intervals produces
#'   different probabilities. For regular data, dividing by a constant time
#'   lag simply rescales uniformly and does not change relative probabilities.
#'   Set to \code{FALSE} only if the data has no meaningful time information.
#'   Zero time lags (duplicate timestamps) are not permitted and will raise
#'   an error -- clean duplicates before running outlier detection.
#' @param silent Logical.  If \code{FALSE} (default) the function
#'   prints a brief running narration (per-iteration counts, threshold
#'   diagnostics, final summary).  Set \code{TRUE} to suppress.  Errors
#'   and warnings are always shown.
#'
#' @return A \code{move2} object. If \code{remove = FALSE}, the following
#'   columns are added:
#'   The added columns are: \code{step_turn_prob} (probability from the 2D
#'   step/turn histogram), \code{delta_step_prob} (probability of the change
#'   in step length), \code{delta_turn_prob} (probability of the change in
#'   turning angle), \code{joint_prob} (product of all three),
#'   \code{outlier_percentile} (0--100, higher = more unusual),
#'   \code{loglr_prob} (signed log-likelihood-ratio of outlier vs not, in
#'   nat units; \code{> 0} where flagged, \code{< 0} below the boundary,
#'   \code{NA} where the detector abstains -- see
#'   \code{DESIGN_evidence_accumulation.md}),
#'   \code{is_outlier} (logical flag), \code{flagged_by_prob} (logical,
#'   same as \code{is_outlier} for this primitive -- named for parity with
#'   the other primitives so the output can be voted by
#'   \code{\link{mt_flag_consensus}}), and \code{is_na_prob} (logical, TRUE
#'   where probability could not be calculated).
#'   If \code{remove = TRUE}, outlier rows (and optionally NA rows) are removed.
#'
#' @examples
#' \dontrun{
#' library(move2)
#'
#' ## load example data
#' fishers <- mt_read(mt_example())
#' fishers <- fishers[!sf::st_is_empty(fishers), ]
#'
#' ## flag outliers per individual (automatic when multiple IDs present)
#' result <- mt_flag_outliers(fishers)
#'
#' ## remove outliers, single individual
#' leroy <- fishers[mt_track_id(fishers) == "M4", ]
#' cleaned <- mt_flag_outliers(leroy, remove = TRUE)
#'
#' ## copula method (faster, works well with small samples)
#' result_cop <- mt_flag_outliers(leroy, method = "copula")
#'
#' ## pooled reference: score against distribution from all individuals
#' result_pop <- mt_flag_outliers(fishers, reference = fishers)
#'
#' ## default: ACF-derived alpha (adapts to the data)
#' result_acf <- mt_flag_outliers(leroy)
#'
#' ## manual override: ignore auto-differences
#' result_a0 <- mt_flag_outliers(leroy, autodiff_alpha = 0)
#'
#' ## manual override: full product (no downweighting)
#' result_a1 <- mt_flag_outliers(leroy, autodiff_alpha = 1.0)
#'
#' ## iterative refinement for consecutive outliers
#' result_iter <- mt_flag_outliers(leroy, iterations = 3)
#'
#' ## quality weighting with satellite count and HDOP
#' result_qc <- mt_flag_outliers(leroy, quality_columns = list(
#'   "gps.satellite.count" = function(s) pnorm(s, mean = 7, sd = 2),
#'   "gps.hdop" = function(h) 1 - pnorm(h, mean = 3, sd = 1.5)
#' ))
#'
#' ## raw step/turn metrics (without time normalisation)
#' result_raw <- mt_flag_outliers(leroy, time_normalize = FALSE)
#' }
#'
#' @references
#' Safi, K. (in preparation). Self-thresholding hierarchical
#' outlier-detection for animal movement tracks. Companion paper to
#' the \pkg{move2utils} R package. Preprint: bioRxiv (DOI forthcoming).
#'
#' @importFrom move2 mt_distance mt_turnangle mt_track_id mt_time mt_time_lags
#' @importFrom sf st_coordinates st_is_empty
#' @importFrom terra rast rasterize ext values merge resample extract vect
#' @importFrom stats ecdf density.default approxfun quantile dweibull optimize pnorm median mad
#' @importFrom graphics par plot points lines legend title
#' @importFrom grDevices colorRampPalette nclass.FD
#' @export
mt_flag_outliers <- function(x, threshold = NULL, prob_type = "joint",
                             remove = FALSE, plot = TRUE, drop_na = FALSE,
                             autodiff_alpha = "acf", method = "histogram",
                             iterations = 1,
                             quality_columns = NULL,
                             time_normalize = TRUE,
                             threshold_type = "gap",
                             step_transform = c("none", "log"),
                             step_floor = 0,
                             reference = NULL,
                             pool_by = NULL,
                             silent = FALSE) {
  step_transform <- match.arg(step_transform)
  ## brief narrator helper -- suppressed under silent = TRUE
  say <- .say(silent)

  ## ---- input validation ----
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  if (!is.numeric(step_floor) || length(step_floor) != 1 ||
      is.na(step_floor) || step_floor < 0) {
    rlang::abort("`step_floor` must be a non-negative scalar (metres).",
                 class = "move2utils_mt_flag_outliers_bad_step_floor")
  }
  ## threshold range validated after default-setting below
  valid_types <- c("joint", "step_turn", "delta_step", "delta_turn", "custom")
  if (!prob_type %in% valid_types) {
    rlang::abort(
      sprintf("`prob_type` must be one of: %s",
              paste(valid_types, collapse = ", ")),
      class = "move2utils_mt_flag_outliers_bad_prob_type")
  }
  ## autodiff_alpha: accept numeric, "acf" (default), or "auto" (legacy)
  acf_alpha <- FALSE
  auto_alpha <- FALSE
  if (is.character(autodiff_alpha) && length(autodiff_alpha) == 1) {
    if (autodiff_alpha == "acf") {
      acf_alpha <- TRUE
    } else if (autodiff_alpha == "auto") {
      auto_alpha <- TRUE
    } else {
      rlang::abort(
        "`autodiff_alpha` must be a non-negative number, \"acf\", or \"auto\".",
        class = "move2utils_mt_flag_outliers_bad_autodiff_alpha")
    }
  } else if (!is.numeric(autodiff_alpha) || length(autodiff_alpha) != 1 ||
             autodiff_alpha < 0) {
    rlang::abort(
      "`autodiff_alpha` must be a non-negative number, \"acf\", or \"auto\".",
      class = "move2utils_mt_flag_outliers_bad_autodiff_alpha")
  }
  valid_methods <- c("histogram", "copula")
  if (!method %in% valid_methods) {
    rlang::abort(
      sprintf("`method` must be one of: %s",
              paste(valid_methods, collapse = ", ")),
      class = "move2utils_mt_flag_outliers_bad_method")
  }
  if (method == "copula") {
    if (!requireNamespace("circular", quietly = TRUE)) {
      rlang::abort(paste0(
        "Package 'circular' is required for method = 'copula'. ",
        "Install it with install.packages('circular')."),
        class = "move2utils_mt_flag_outliers_missing_circular")
    }
    if (!requireNamespace("MASS", quietly = TRUE)) {
      rlang::abort(paste0(
        "Package 'MASS' is required for method = 'copula'. ",
        "Install it with install.packages('MASS')."),
        class = "move2utils_mt_flag_outliers_missing_MASS")
    }
  }
  if (!is.numeric(iterations) || length(iterations) != 1 || iterations < 1) {
    rlang::abort("`iterations` must be a positive integer.",
                 class = "move2utils_mt_flag_outliers_bad_iterations")
  }
  iterations <- as.integer(iterations)
  if (!is.null(quality_columns)) {
    if (!is.list(quality_columns) || is.null(names(quality_columns))) {
      rlang::abort("`quality_columns` must be a named list of functions.",
                   class = "move2utils_mt_flag_outliers_bad_quality_columns")
    }
    for (nm in names(quality_columns)) {
      if (!nm %in% names(x)) {
        rlang::abort(
          sprintf("Column '%s' not found in `x`.", nm),
          class = "move2utils_mt_flag_outliers_missing_quality_column")
      }
      if (!is.function(quality_columns[[nm]])) {
        rlang::abort(
          sprintf("Element '%s' in `quality_columns` must be a function.", nm),
          class = "move2utils_mt_flag_outliers_bad_quality_column_fn")
      }
    }
  }
  valid_threshold_types <- c("gap", "significance", "percentile", "entropy")
  if (!threshold_type %in% valid_threshold_types) {
    rlang::abort(
      sprintf("`threshold_type` must be one of: %s",
              paste(valid_threshold_types, collapse = ", ")),
      class = "move2utils_mt_flag_outliers_bad_threshold_type")
  }
  ## set default threshold based on threshold_type
  ##
  ## Each default below is HEURISTIC -- defensible but not derived from
  ## first principles. Sensible alternative ranges noted alongside.
  ##
  ##   gap = 3
  ##     Break-size multiplier vs. local noise in the broken-stick
  ##     score. The 3-sigma convention applied loosely. Plausible
  ##     range: 2--5.  Tighter values increase false-positive rate on
  ##     unimodal tails; looser values miss subtle outliers.
  ##
  ##   significance = 0.001 / percentile = 0.001
  ##     Bottom 0.1% of the joint-probability distribution. Plausible
  ##     range: 1e-4 -- 1e-2. Tighter values are more conservative.
  ##
  ##   entropy = 0.3
  ##     Unified package-wide entropy default (2026-05-09). Maximum
  ##     valley-to-peak density ratio for log-probability valleys;
  ##     same numerical value used by the bridge and speed-cap
  ##     primitives. Validated by the 2026-05-06 Raven sensitivity
  ##     sweep on the cascade-relevant primitives; the prob
  ##     primitive's entropy mode is not exercised by the cascade
  ##     (cascade default is gap), so the unification touches only
  ##     standalone use. Plausible range: 0.3--0.7.
  ## NULL threshold defers to the leaf detectors' formal defaults
  ## (single source of truth, after the 2026-05-25 propagation fix).
  ##
  ## - For threshold_type = "gap" we forward NULL through to
  ##   .gap_threshold_lower's formal at L763.
  ## - For threshold_type = "entropy" the entropy path here is an
  ##   inline KDE-valley implementation (not a call to
  ##   .entropy_threshold_lower) so we resolve NULL by reading the
  ##   leaf's formal default at runtime (sweep-override-friendly).
  ## - For "significance" and "percentile" the threshold is consumed
  ##   inline; resolve to the established 0.001 default.
  if (is.null(threshold)) {
    threshold <- switch(threshold_type,
      gap          = NULL,   # forward NULL; leaf formal is the truth
      entropy      = formals(.entropy_threshold_lower)$threshold,
      significance = 0.001,
      percentile   = 0.001
    )
  }
  if (!is.null(threshold) &&
      (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0)) {
    rlang::abort("`threshold` must be a positive number or NULL.",
                 class = "move2utils_mt_flag_outliers_bad_threshold")
  }
  if (!is.logical(time_normalize) || length(time_normalize) != 1) {
    rlang::abort("`time_normalize` must be TRUE or FALSE.",
                 class = "move2utils_mt_flag_outliers_bad_time_normalize")
  }
  if (!is.null(reference) && !inherits(reference, "move2")) {
    rlang::abort("`reference` must be a move2 object or NULL.",
                 class = "move2utils_mt_flag_outliers_bad_reference")
  }
  if (!is.null(reference) && !is.null(pool_by)) {
    rlang::abort(paste0(
      "`reference =` and `pool_by =` are mutually exclusive. ",
      "`reference` fits a single global distribution from an external clean track. ",
      "`pool_by` fits one distribution per group within the input. ",
      "Choose one."),
      class = "move2utils_mt_flag_outliers_reference_with_pool_by")
  }
  if (!is.null(pool_by)) {
    ## Shape, column existence, and (for length-2 input) strict
    ## nesting are validated by `.resolve_pool_groups`.
    invisible(.resolve_pool_groups(x, pool_by, silent = TRUE))
  }

  ## ---- iterative refinement ----
  ## Iteration > 1 takes the legacy path which recursively calls
  ## mt_flag_outliers().  Within the loop, iterations is fixed to 1 so
  ## each inner call goes through the dispatcher-based fast path.
  if (iterations > 1) {
    return(.iterative_flag(x, threshold = threshold, prob_type = prob_type,
                           remove = remove, plot = plot, drop_na = drop_na,
                           autodiff_alpha = if (acf_alpha) "acf" else if (auto_alpha) "auto" else autodiff_alpha,
                           method = method,
                           iterations = iterations,
                           quality_columns = quality_columns,
                           time_normalize = time_normalize,
                           threshold_type = threshold_type,
                           step_transform = step_transform,
                           step_floor = step_floor,
                           reference = reference,
                           pool_by = pool_by,
                           silent = silent))
  }

  ## ---- precompute reference distributions if external reference ----
  ## In cascade scope and standalone scope without `reference` and
  ## without `pool_by`, the reference is the input itself -- .fn_core
  ## handles that case locally.  Only the external-reference case
  ## needs precomputation here (single global ref), since the
  ## dispatcher does not see `reference`.
  ref_args <- NULL
  if (!is.null(reference)) {
    say("Building distributions from reference data...")
    ref_args <- .build_prob_ref_args(reference, time_normalize)
  }

  ## ---- precompute pool_by per-group references --------------------
  ## When pool_by is set, fit ref_step / ref_turn / ref_delta_* once
  ## per pool group from the union of that group's tracks, then build
  ## a `pool_args_by_track` map injecting the group's fitted args
  ## into each track's dispatcher call.
  pool_args_by_track <- NULL
  if (!is.null(pool_by)) {
    maps <- .resolve_pool_groups(x, pool_by, silent = silent)
    ## mt_flag_outliers uses an *integrated* pool: the ref distribution
    ## is fit from the OUTER group and injected per-track.  The INNER
    ## map has no role for this primitive (there is no post-hoc union
    ## step where it would be the operating unit) -- under length-2
    ## `pool_by`, inner is harmlessly unused.  Validation already ran
    ## inside `.resolve_pool_groups` (nesting check, column existence,
    ## etc.) so the inner-as-unused path is still safe.
    outer_map <- maps$outer
    n_groups <- length(unique(outer_map))
    say(sprintf(
      "Pooling by [%s] -- %d outer group(s) across %d track(s).",
      paste(pool_by, collapse = ","), n_groups, length(outer_map)))
    ids_event <- as.character(move2::mt_track_id(x))
    group_fits <- list()
    for (g in unique(outer_map)) {
      g_tracks <- names(outer_map)[outer_map == g]
      g_idx    <- which(ids_event %in% g_tracks)
      x_g      <- x[g_idx, , drop = FALSE]
      group_fits[[g]] <- .build_prob_ref_args(x_g, time_normalize)
    }
    pool_args_by_track <- lapply(outer_map, function(g) group_fits[[g]])
    names(pool_args_by_track) <- names(outer_map)
  }

  ## ---- precompute quality weight if quality_columns supplied ----
  quality_weight <- NULL
  if (!is.null(quality_columns)) {
    n <- nrow(x)
    quality_weight <- rep(1.0, n)
    for (nm in names(quality_columns)) {
      raw_vals <- as.numeric(x[[nm]])
      qfun <- quality_columns[[nm]]
      w <- as.numeric(qfun(raw_vals))
      w <- pmin(pmax(w, 0), 1)
      w[is.na(w)] <- 1.0
      quality_weight <- quality_weight * w
    }
  }

  was_longlat <- isTRUE(sf::st_is_longlat(x))

  ## ---- common plumbing: multi-track dispatch + extraction + hygiene
  ## + .prob_fn_core + lift.  Stays in input CRS (Haversine + spherical
  ## azimuth on lon/lat, Euclidean on projected) -- the prob math is
  ## CRS-invariant within sub-pct accuracy via the helpers. ----
  out <- .clean_track_dispatch(
    x,
    fn_core             = .prob_fn_core,
    fn_core_args        = c(
      list(was_longlat    = was_longlat,
            threshold      = threshold,
            prob_type      = prob_type,
            autodiff_alpha = autodiff_alpha,
            acf_alpha      = acf_alpha,
            auto_alpha     = auto_alpha,
            method         = method,
            time_normalize = time_normalize,
            threshold_type = threshold_type,
            step_transform = step_transform,
            step_floor     = step_floor,
            silent         = silent),
      ref_args),
    per_track_args      = list(quality_weight = quality_weight),
    pool_args_by_track  = pool_args_by_track,
    need_time           = TRUE,
    hygiene_strict_time = FALSE,
    ## canonicalise to AEQD: the probability surface uses turn angles, which
    ## are CRS-dependent (geographic vs Cartesian), so the input CRS must not
    ## leak into the result (was FALSE, which made longlat != projected).
    project_longlat     = TRUE,
    n_min               = 3L,
    n_min_severity      = "say",
    primitive_label     = "prob",
    silent              = silent)

  ## ---- diagnostic plot ----
  if (plot) {
    prob_col <- switch(prob_type,
      joint = "joint_prob", step_turn = "step_turn_prob",
      delta_step = "delta_step_prob", delta_turn = "delta_turn_prob",
      custom = "custom_prob")
    prob_label <- switch(prob_type,
      joint = "Joint probability", step_turn = "Step-turn probability",
      delta_step = "Delta step probability",
      delta_turn = "Delta turn probability",
      custom = "Custom probability (step-turn x delta-step)")
    .plot_outliers(out, prob_col, prob_label, out$is_outlier,
                   out$is_na_prob, drop_na,
                   sum(out$is_outlier, na.rm = TRUE),
                   sum(out$is_na_prob))
  }

  ## ---- log-LR emission (evidence currency; additive, is_outlier untouched) ----
  out <- .attach_prob_loglr(out)

  ## ---- return ----
  if (remove) {
    to_remove <- out$is_outlier
    if (drop_na) to_remove <- to_remove | out$is_na_prob
    return(out[!to_remove, ])
  }
  out
}


## Raw-matrix entry point used by mt_clean_track to skip per-iter
## sf-class slicing.  Faithful port of the single-iteration prob math
## that used to live in `mt_flag_outliers()`'s wrapper.  Operates on
## raw matrices and indexed-by-active_idx; returns active-indexed
## result vectors.
##
## Quality columns are pre-resolved to a length-n `quality_weight`
## vector by the wrapper (the .fn_core knows nothing about move2
## columns).  External reference distributions, when supplied, are
## precomputed by the wrapper and passed as `ref_*` args.
##
## @keywords internal
.prob_fn_core <- function(cc, t_s, active_idx, was_longlat,
                            threshold, prob_type, autodiff_alpha,
                            acf_alpha, auto_alpha,
                            method, time_normalize, threshold_type,
                            step_transform, step_floor,
                            ref_step       = NULL,
                            ref_turn       = NULL,
                            ref_deltaStep  = NULL,
                            ref_deltaTurn  = NULL,
                            ref_delta_gaps = NULL,
                            quality_weight = NULL,
                            silent         = TRUE) {

  ## Contract: active_idx sorted ascending with no duplicates.
  stopifnot(!is.unsorted(active_idx), !anyDuplicated(active_idx))

  say <- function(...) if (!silent) message(...)
  n_a <- length(active_idx)
  cc_a <- cc[active_idx, , drop = FALSE]
  t_a  <- t_s[active_idx]

  say("Calculating movement metrics...")
  stepLength <- .step_lengths_from_cc(cc_a, was_longlat)
  turnAngle  <- .turn_angles_from_cc (cc_a, was_longlat)
  turnAngle_raw <- turnAngle

  if (sum(!is.na(stepLength) & !is.na(turnAngle)) < 3) {
    rlang::abort(
      "Not enough valid locations to compute probabilities (need >= 3).",
      class = "move2utils_mt_flag_outliers_too_few_valid_locations")
  }

  ## ---- time normalisation ----
  if (time_normalize) {
    time_lags <- c(diff(t_a), NA_real_)
    n_zero <- sum(time_lags == 0, na.rm = TRUE)
    if (n_zero > 0) {
      rlang::abort(
        sprintf(paste0(
          "Found %d zero time lag(s) (duplicate timestamps). ",
          "Remove duplicates before running outlier detection -- see ",
          "move2::mt_filter_unique() or the cleaning steps in chapter 2."),
          n_zero),
        class = c("move2utils_input_duplicate_timestamps",
                  "move2utils_mt_flag_outliers_zero_time_lag"))
    }
    na_tl <- is.na(time_lags)
    speed            <- ifelse(na_tl, NA_real_, stepLength / time_lags)
    angular_velocity <- ifelse(na_tl, NA_real_, turnAngle  / time_lags)
    stepLength <- speed
    turnAngle  <- angular_velocity
  } else {
    time_lags <- NULL
  }

  deltaStep <- c(NA, diff(stepLength))
  deltaTurn <- c(NA, diff(turnAngle))

  ## ---- ACF-derived alpha ----
  if (acf_alpha) {
    spd_valid <- stepLength[!is.na(stepLength)]
    ang_valid <- turnAngle [!is.na(turnAngle)]
    if (length(spd_valid) >= 10 && length(ang_valid) >= 10) {
      r_v <- max(stats::acf(spd_valid, lag.max = 1, plot = FALSE)$acf[2], 0)
      r_w <- max(stats::acf(ang_valid, lag.max = 1, plot = FALSE)$acf[2], 0)
      autodiff_alpha <- sqrt(r_v * r_w)
      say(sprintf(
        "ACF-derived alpha: %.3f (r_speed=%.3f, r_angvel=%.3f)",
        autodiff_alpha, r_v, r_w))
    } else {
      autodiff_alpha <- 0.5
      say("Too few valid values for ACF estimation; using alpha = 0.5")
    }
  }

  ## ---- gap lengths for autodifference scaling ----
  if (time_normalize) {
    delta_gaps <- c(NA, (time_lags[-length(time_lags)] + time_lags[-1]) / 2)
  } else {
    delta_gaps <- rep(1, length(stepLength))
  }

  ## ---- reference distributions ----
  if (is.null(ref_step)) {
    ref_step       <- stepLength
    ref_turn       <- turnAngle
    ref_deltaStep  <- deltaStep
    ref_deltaTurn  <- deltaTurn
    ref_delta_gaps <- delta_gaps
  } else {
    say("Building distributions from supplied reference data...")
  }

  say("Calculating probability distributions...")
  step_turn_prob  <- rep(NA_real_, n_a)
  delta_step_prob <- rep(NA_real_, n_a)
  delta_turn_prob <- rep(NA_real_, n_a)
  joint_prob      <- rep(NA_real_, n_a)

  ga_step <- .gap_aware_autodiff(ref_deltaStep, ref_delta_gaps)
  ga_turn <- .gap_aware_autodiff(ref_deltaTurn, ref_delta_gaps)

  ## Surface whether the gap-aware auto-difference is actually in force.
  ## When it falls back to a constant scale the numerics stay valid but
  ## the detector is no longer gap-aware, and nothing downstream can tell.
  for (nm in c(step = "step", turn = "turning-angle")) {
    ga <- if (nm == "step") ga_step else ga_turn
    if (isFALSE(ga$gap_dependent)) {
      say(sprintf(
        "Note: gap-aware scaling is NOT in force for the %s auto-difference (%s);\n  the scale is constant, so autodifferences are not gap-normalised.",
        nm, ga$gap_reason))
    }
  }

  if (method == "histogram") {
    if (step_transform == "none") {
      valid_step <- ref_step[!is.na(ref_step) & ref_step > 0]
      if (length(valid_step) >= 10) {
        iqr <- stats::IQR(valid_step)
        rng <- diff(range(valid_step))
        if (iqr > 0 && rng / iqr > 500) {
          say(sprintf(
            "Note: step-length range/IQR = %.0f is extreme;\n",
            rng / iqr),
            "  teleport-class GPS errors are better handled by\n",
            "  mt_filter_gps_quality() (drop fixes with <5 satellites)\n",
            "  and mt_flag_outliers_bridge() (geometric, leverage-immune).\n",
            "  step_transform = \"log\" is available but can hide\n",
            "  physiologically-plausible joint turn/step outliers."
          )
        }
      }
    }
    step_fn <- switch(step_transform, log = log1p, none = identity)
    ref_step_h <- step_fn(ref_step)
    step_h     <- step_fn(stepLength)

    hist2d <- .turn_step_hist(ref_turn, ref_step_h)

    if (is.null(hist2d)) {
      ## Degenerate reference: no positive step lengths (e.g. an
      ## all-stationary state subset of repeated positions).  No
      ## (turn, step) histogram can be built; leave step_turn_prob as
      ## NA so these fixes are kept rather than flagged as kinematic
      ## outliers.
      step_turn_prob <- rep(NA_real_, n_a)
    } else {
      coords_st <- cbind(turnAngle, step_h)
      extracted <- terra::extract(hist2d, coords_st)
      step_turn_prob <- extracted[, 1]
      ok_stp <- !is.na(step_turn_prob)
      step_turn_prob[ok_stp] <- pmax(step_turn_prob[ok_stp], .Machine$double.xmin)
    }
  } else {
    step_turn_prob <- .parametric_step_turn_prob(stepLength, turnAngle_raw)
  }

  delta_step_prob <- .gap_aware_delta_prob(deltaStep, delta_gaps, ga_step,
                                           wrap = FALSE)
  delta_turn_prob <- .gap_aware_delta_prob(deltaTurn, delta_gaps, ga_turn,
                                           wrap = TRUE)

  if (auto_alpha) {
    autodiff_alpha <- .optimize_alpha(step_turn_prob, delta_step_prob,
                                      delta_turn_prob)
    say(sprintf("Auto-optimised alpha: %.4f", autodiff_alpha))
  }

  say("Calculating joint probabilities...")
  valid <- !is.na(step_turn_prob) & !is.na(delta_step_prob) &
           !is.na(delta_turn_prob)
  joint_prob[valid] <- step_turn_prob[valid] *
    (delta_step_prob[valid] * delta_turn_prob[valid])^autodiff_alpha

  ## ---- apply quality weighting ----
  qw_active <- NULL
  if (!is.null(quality_weight)) {
    qw_active <- quality_weight[active_idx]
    joint_prob <- joint_prob * qw_active
  }

  custom_prob <- if (prob_type == "custom") {
    step_turn_prob * delta_step_prob
  } else NULL

  probs <- switch(prob_type,
    joint      = joint_prob,
    step_turn  = step_turn_prob,
    delta_step = delta_step_prob,
    delta_turn = delta_turn_prob,
    custom     = custom_prob)
  prob_label <- switch(prob_type,
    joint      = "Joint probability",
    step_turn  = "Step-turn probability",
    delta_step = "Delta step probability",
    delta_turn = "Delta turn probability",
    custom     = "Custom probability (step-turn x delta-step)")

  ## ---- identify outliers ----
  say("Identifying outliers...")
  is_na_prob <- is.na(probs)
  outlier_percentile <- rep(NA_real_, n_a)
  outlier_zscore     <- rep(NA_real_, n_a)
  is_outlier         <- rep(FALSE, n_a)

  if (threshold_type == "gap") {
    probs_floored <- pmax(probs, .Machine$double.xmin)
    log_probs_all <- log(probs_floored)
    log_probs_all[is_na_prob] <- NA_real_
    log_probs_valid <- log_probs_all[!is_na_prob & is.finite(log_probs_all)]

    if (length(log_probs_valid) < 10) {
      rlang::warn(paste0(
        "Too few valid probabilities for gap-based detection. ",
        "Falling back to percentile method."),
        class = "move2utils_mt_flag_outliers_too_few_for_gap")
      threshold_type <- "percentile"
    } else {
      ## NULL threshold defers to .gap_threshold_lower's leaf formal
      ## (single source of truth).  See
      ## audits/2026-05-25-parameter-propagation/findings.md §1.
      gap_res <- if (is.null(threshold))
                   .gap_threshold_lower(log_probs_all)
                 else
                   .gap_threshold_lower(log_probs_all, threshold = threshold)
      is_outlier <- gap_res$is_outlier & !is_na_prob
      is_outlier[is.na(is_outlier)] <- FALSE
      outlier_percentile <- gap_res$percentile
    }
  }

  if (threshold_type == "significance") {
    log_probs <- log(pmax(probs[!is_na_prob], .Machine$double.xmin))
    log_probs <- log_probs[is.finite(log_probs)]

    if (length(log_probs) < 5) {
      rlang::warn(paste0(
        "Too few valid probabilities for significance-based thresholding. ",
        "Falling back to percentile method."),
        class = "move2utils_mt_flag_outliers_too_few_for_significance")
      threshold_type <- "percentile"
    } else {
      bulk_median <- stats::median(log_probs)
      bulk_mad    <- stats::mad(log_probs)

      if (bulk_mad == 0) {
        outlier_zscore     <- rep(0, n_a)
        outlier_percentile <- rep(50, n_a)
        is_outlier         <- rep(FALSE, n_a)
      } else {
        z <- (log(pmax(probs, .Machine$double.xmin)) - bulk_median) / bulk_mad
        p_val <- stats::pnorm(z)
        is_outlier <- p_val < threshold & !is_na_prob
        is_outlier[is.na(is_outlier)] <- FALSE
        outlier_zscore     <- round(z, 4)
        outlier_percentile <- round(p_val * 100, 4)
      }
    }
  }

  if (threshold_type == "entropy") {
    probs_floored <- pmax(probs, .Machine$double.xmin)
    log_probs_all <- log(probs_floored)
    log_probs_all[is_na_prob] <- NA_real_
    log_probs_valid <- log_probs_all[!is_na_prob & is.finite(log_probs_all)]

    if (length(log_probs_valid) < 10) {
      rlang::warn(paste0(
        "Too few valid probabilities for entropy-based detection. ",
        "Falling back to percentile method."),
        class = "move2utils_mt_flag_outliers_too_few_for_entropy")
      threshold_type <- "percentile"
    } else {
      d <- stats::density(log_probs_valid, n = 512)
      peak_idx <- which.max(d$y)
      peak_density <- d$y[peak_idx]

      lower <- which(d$x < d$x[peak_idx])
      valley_x <- NULL
      valley_depth <- Inf

      if (length(lower) >= 3) {
        for (i in 2:(length(lower) - 1)) {
          j <- lower[i]
          if (d$y[j] < d$y[j - 1] && d$y[j] < d$y[j + 1]) {
            depth_ratio <- d$y[j] / peak_density
            if (depth_ratio < valley_depth) {
              valley_depth <- depth_ratio
              valley_x <- d$x[j]
            }
          }
        }
      }

      if (!is.null(valley_x) && valley_depth < threshold) {
        is_outlier <- log_probs_all < valley_x & !is_na_prob
        is_outlier[is.na(is_outlier)] <- FALSE
      } else {
        is_outlier <- rep(FALSE, n_a)
      }

      outlier_percentile <- round(
        100 * stats::ecdf(log_probs_valid)(log_probs_all), 4
      )
    }
  }

  if (threshold_type == "percentile") {
    ## When we ARRIVE here as a fallback from gap/entropy/significance
    ## (too few valid probabilities for that method), `threshold` still
    ## carries that method's value, which is NOT a percentile fraction:
    ## it may be NULL (deferred to a leaf formal -- `100 - NULL * 100`
    ## is numeric(0), giving a zero-length `is_outlier` that then breaks
    ## the caller's `x[idx] <- is_outlier` scatter), or a gap-spacing
    ## multiple like 3 (giving a nonsensical `100 - 300`).  Coerce
    ## anything outside the valid percentile-fraction range (0, 1] to the
    ## percentile default.
    pct_threshold <- if (is.null(threshold) || length(threshold) != 1L ||
                         !is.finite(threshold) || threshold <= 0 || threshold > 1) {
      0.001
    } else threshold
    valid_probs <- probs[!is_na_prob]
    if (length(valid_probs) == 0L) {
      ## No valid probabilities at all (e.g. an all-stationary subset
      ## whose step-turn histogram was degenerate, leaving every
      ## joint_prob NA).  stats::ecdf() errors on empty input; there is
      ## nothing to rank, so flag nothing.
      is_outlier <- rep(FALSE, n_a)
    } else {
      ecdf_func <- stats::ecdf(valid_probs)
      outlier_percentile <- round((1 - ecdf_func(probs)) * 100, 2)
      is_outlier <- outlier_percentile >= (100 - pct_threshold * 100)
      is_outlier[is.na(is_outlier)] <- FALSE
    }
  }

  if (step_floor > 0) {
    n_before <- sum(is_outlier, na.rm = TRUE)
    under_floor <- !is.na(stepLength) & stepLength <= step_floor
    is_outlier <- is_outlier & !under_floor
    n_downgraded <- n_before - sum(is_outlier, na.rm = TRUE)
    if (n_downgraded > 0) {
      say(sprintf(
        "step_floor = %g m: downgraded %d rate-flagged fixes whose absolute step was within GPS noise.",
        step_floor, n_downgraded))
    }
  }

  ## ---- report ----
  n_outliers <- sum(is_outlier, na.rm = TRUE)
  n_valid    <- sum(!is_na_prob)
  n_na       <- sum(is_na_prob)

  if (threshold_type == "gap") {
    say(sprintf(
      "Identified %d outliers (%.2f%% of %d locations) based on %s.\nGap threshold: %.1f x local spacing (non-parametric break detection).",
      n_outliers, 100 * n_outliers / n_valid, n_valid, prob_label, threshold
    ))
  } else if (threshold_type == "significance") {
    say(sprintf(
      "Identified %d outliers (%.2f%% of %d locations) based on %s.\nSignificance threshold: %.4g (robust z-score on log-probabilities).",
      n_outliers, 100 * n_outliers / n_valid, n_valid, prob_label, threshold
    ))
  } else if (threshold_type == "entropy") {
    say(sprintf(
      "Identified %d outliers (%.2f%% of %d locations) based on %s.\nEntropy threshold: density-valley detection (max valley/peak ratio: %.2f).",
      n_outliers, 100 * n_outliers / n_valid, n_valid, prob_label, threshold
    ))
  } else {
    say(sprintf(
      "Identified %d outliers (%.2f%% of %d locations) based on %s.\nPercentile threshold: %s (%.4g).",
      n_outliers, 100 * n_outliers / n_valid, n_valid, prob_label,
      paste0(threshold * 100, "th"), threshold
    ))
  }
  if (n_na > 0) {
    n_stationary <- sum(stepLength == 0, na.rm = TRUE)
    na_reason <- if (n_stationary > 0) {
      sprintf(" (includes %d stationary fixes)", n_stationary)
    } else ""
    say(sprintf(
      "%d locations (%.1f%%) have NA probabilities%s --will be kept.",
      n_na, 100 * n_na / n_a, na_reason))
  }

  out_list <- list(
    is_outlier         = is_outlier,
    is_na_prob         = is_na_prob,
    step_turn_prob     = step_turn_prob,
    delta_step_prob    = delta_step_prob,
    delta_turn_prob    = delta_turn_prob,
    joint_prob         = joint_prob,
    outlier_percentile = outlier_percentile
  )
  if (prob_type == "custom") {
    out_list$custom_prob <- custom_prob
  }
  if (threshold_type == "significance") {
    out_list$outlier_zscore <- outlier_zscore
  }
  if (!is.null(qw_active)) {
    out_list$quality_weight <- qw_active
  }
  out_list
}


# ---- internal helpers (not exported) ----

#' Wrap angle to the range from -pi to pi
#' @noRd
.wrap_angle <- function(x) {
  (x + pi) %% (2 * pi) - pi
}

#' Gap-normalised delta probability, vectorised over the whole track
#'
#' Drop-in replacement for the per-row loop that used to call
#' \code{scale_fun} and \code{kde_fun} once per location. Both underlying
#' functions are \code{approxfun}-based and vector-safe, so this runs in
#' a single pass. Index 1 is always \code{NA} (no prior step), matching
#' the old loop that started at \code{i = 2}.
#'
#' @param delta Numeric vector of autodifferences (length n).
#' @param gaps Numeric vector of gap lengths (length n).
#' @param ga The \code{.gap_aware_autodiff()} result (list with
#'   \code{scale_fun}, \code{kde_fun}).
#' @param wrap Logical. If \code{TRUE}, wrap \code{delta} into
#'   \eqn{[-\pi, \pi)} before evaluating (used for the angular branch).
#' @return Numeric vector of length \code{n}.
#' @noRd
.gap_aware_delta_prob <- function(delta, gaps, ga, wrap = FALSE) {
  n <- length(delta)
  out <- rep(NA_real_, n)
  if (n < 2) return(out)

  ok <- !is.na(delta) & !is.na(gaps)
  ok[1] <- FALSE
  if (!any(ok)) return(out)

  s <- ga$scale_fun(gaps[ok])
  pos <- !is.na(s) & s > 0
  if (!any(pos)) return(out)

  idx <- which(ok)[pos]
  s_pos <- s[pos]
  d <- delta[idx]
  if (wrap) d <- .wrap_angle(d)
  out[idx] <- pmax(ga$kde_fun(d / s_pos) / s_pos, .Machine$double.xmin)
  out
}

#' Safe kernel density interpolation function
#'
#' Wraps density() + approxfun() with guards for small samples and
#' protection against negative extrapolation. Returns a function that
#' always produces non-negative values (floored at .Machine$double.xmin).
#' @param vals Numeric vector (NAs already removed by caller).
#' @return A function that maps numeric values to density estimates.
#' @noRd
.safe_kde_fun <- function(vals) {
  if (length(vals) < 3) {
    ## too few values for meaningful KDE — return uniform density
    rng <- range(vals)
    width <- max(rng[2] - rng[1], 1)
    return(function(x) rep(1 / width, length(x)))
  }
  kde <- stats::density.default(vals)
  fn <- stats::approxfun(kde, rule = 2)  ## rule=2: constant extrapolation
  ## wrap to guarantee non-negative output
  function(x) pmax(fn(x), .Machine$double.xmin)
}

#' Estimate gap-dependent scale of autodifferences
#'
#' For irregularly sampled data, the expected magnitude of speed (or
#' angular velocity) changes depends on the time gap between steps.
#' This function estimates the relationship non-parametrically from the
#' data, returning a function that maps gap length to the expected
#' scale (MAD) of the autodifference.
#'
#' The estimate is built by binning (|delta|, gap) pairs into
#' equal-count bins by gap, computing the MAD in each bin, and
#' interpolating.  This adapts to any movement process without
#' assuming a specific model.
#'
#' @param deltas Numeric vector of autodifferences (Δv or Δω).
#' @param gaps Numeric vector of gap lengths (same units as time lags).
#'   Must be the same length as \code{deltas}.
#' @param n_bins Number of quantile-based bins.  Default 8 is
#'   HEURISTIC -- a round number balancing local-scale resolution
#'   (more bins = finer gap-dependent MAD estimate) against
#'   per-bin sample size.  Plausible range: 5--15.  The function
#'   additionally clamps n_bins so each bin holds at least ~5
#'   fixes, with a minimum of 2 bins.
#' @return A list with four elements:
#'   \item{scale_fun}{A function mapping gap → expected scale (MAD).}
#'   \item{kde_fun}{A KDE function on the gap-normalised deltas
#'     (delta / scale).  Returns density in the normalised space.}
#'   \item{gap_dependent}{Logical.  \code{TRUE} when the returned
#'     \code{scale_fun} actually varies with the gap; \code{FALSE} when
#'     the estimator fell back to a constant scale.  The fallback is
#'     silent to the numerics but changes what the detector is doing,
#'     so it is reported rather than left for the caller to infer --
#'     the gap-aware auto-difference is the package's headline
#'     contribution and "it is not in force here" is information the
#'     user is entitled to.}
#'   \item{gap_reason}{Character.  Why gap-dependence was or was not
#'     established, suitable for narration.}
#' @noRd
.gap_aware_autodiff <- function(deltas, gaps, n_bins = 8) {
  ok <- !is.na(deltas) & !is.na(gaps) & gaps > 0
  d <- deltas[ok]
  g <- gaps[ok]

  if (length(d) < 10) {
    ## too few points — fall back to gap-unaware KDE
    kde <- .safe_kde_fun(d)
    return(list(
      scale_fun = function(gap) rep(1, length(gap)),
      kde_fun   = kde,
      gap_dependent = FALSE,
      gap_reason = sprintf(
        "only %d valid (delta, gap) pair(s); need 10", length(d))
    ))
  }

  ## quantile-based bins (equal count, robust to skewed gap distribution).
  ## Soft floor of 5 fixes per bin (HEURISTIC; plausible 3--10):
  ## prevents tiny bins from producing degenerate MAD estimates.
  ## Hard floor of 2 bins so the gap-dependence is at least binary.
  n_bins <- min(n_bins, length(d) %/% 5)
  n_bins <- max(n_bins, 2)
  breaks <- stats::quantile(g, probs = seq(0, 1, length.out = n_bins + 1))
  ## ensure unique breaks
  breaks <- unique(breaks)
  if (length(breaks) < 3) {
    ## Fewer than 3 distinct quantile breaks -> the binning cannot resolve
    ## a gap trend, so fall back to a constant scale.  NOTE this is NOT
    ## the same as "all gaps are the same" (which is what this comment
    ## used to claim).  A strictly two-valued gap distribution -- the
    ## ordinary day/night duty cycle -- can land here or on the binned
    ## path depending only on the MIXING PROPORTION of the two values,
    ## because that decides whether any quantile probe interpolates
    ## strictly between them.  So the same tag on the same schedule can
    ## change regime between seasons.  Reported via `gap_dependent`.
    kde <- .safe_kde_fun(d)
    return(list(
      scale_fun = function(gap) rep(1, length(gap)),
      kde_fun   = kde,
      gap_dependent = FALSE,
      gap_reason = sprintf(
        paste0("gap distribution yields only %d distinct quantile ",
               "break(s); need 3 (a strictly two-valued duty cycle can ",
               "land here depending on its mixing proportion)"),
        length(breaks))
    ))
  }
  breaks[1] <- breaks[1] - 1e-10
  bins <- cut(g, breaks = breaks, include.lowest = TRUE)

  ## MAD in each bin
  bin_mid <- tapply(g, bins, stats::median)
  bin_mad <- tapply(d, bins, stats::mad)
  ## replace zero MADs with the overall MAD
  overall_mad <- stats::mad(d)
  if (overall_mad == 0) overall_mad <- stats::sd(d) * 0.6745
  if (overall_mad == 0) overall_mad <- 1
  bin_mad[is.na(bin_mad) | bin_mad == 0] <- overall_mad

  ## interpolation function: gap → expected MAD
  bin_mid_v <- as.numeric(bin_mid)
  bin_mad_v <- as.numeric(bin_mad)
  ord <- order(bin_mid_v)
  scale_fun <- stats::approxfun(bin_mid_v[ord], bin_mad_v[ord],
                                 rule = 2)  # constant extrapolation

  ## normalise deltas and build KDE
  predicted_scale <- scale_fun(g)
  predicted_scale[predicted_scale <= 0] <- overall_mad
  d_norm <- d / predicted_scale

  kde <- .safe_kde_fun(d_norm)

  ## A binned path can still be constant in effect: if every bin's MAD was
  ## zero or NA it was replaced by `overall_mad` above, and approxfun on
  ## identical ordinates is a constant function.  Report that as not
  ## gap-dependent too -- what matters to the caller is whether the scale
  ## varies with the gap, not which branch produced it.
  varies <- length(unique(bin_mad_v)) > 1L
  list(scale_fun = scale_fun, kde_fun = kde,
       gap_dependent = varies,
       gap_reason = if (varies) {
         sprintf("estimated across %d gap bin(s)", length(bin_mad_v))
       } else {
         sprintf(paste0("%d gap bin(s) formed but all bin scales are ",
                        "identical; scale is constant in effect"),
                 length(bin_mad_v))
       })
}


#' Build 2D histogram of turning angles (circular) and step lengths
#'
#' Returns a standardised terra SpatRaster suitable for probability extraction.
#' @noRd
.turn_step_hist <- function(turn_angle, step_length) {
  ## Guard degenerate inputs before the FD bin-count rule (which errors
  ## on an empty vector) and before building the raster (terra rejects a
  ## zero-height ymin == ymax extent).  This happens on an all-stationary
  ## state subset of repeated positions (every step length 0), or a
  ## subset with no valid turn angles -- the (turn, step) histogram is
  ## undefined there.  Signal the caller to skip the histogram path.
  valid_turn <- turn_angle[!is.na(turn_angle)]
  valid_step <- step_length[!is.na(step_length)]
  max_step   <- if (length(valid_step)) max(valid_step) else NA_real_
  if (length(valid_turn) == 0L || !is.finite(max_step) || max_step <= 0) {
    return(NULL)
  }

  ## Determine bin counts using Freedman-Diaconis, but cap at the
  ## smoothing target size. On tracks with a heavy-tailed step
  ## distribution (e.g. a local-foraging bird plus a handful of
  ## teleport-class GPS errors) the FD rule blows up the range / IQR
  ## ratio, producing rasters with O(n^2) cells. The downstream
  ## bilinear resample to 150x150 is lossy anyway, so there is no
  ## information upside to building bigger than that; there is a
  ## large runtime downside.
  max_bins <- 200L
  nx <- max(min(nclass.FD(valid_turn), max_bins), 12L)
  ny <- max(min(nclass.FD(valid_step), max_bins), 12L)

  ymax <- 1.1 * max_step

  ## rasterise the 2D histogram
  r <- terra::rast(ncol = nx, nrow = ny,
                   xmin = -pi, xmax = pi,
                   ymin = 0, ymax = ymax,
                   crs = NA)
  xy <- cbind(as.numeric(turn_angle), as.numeric(step_length))
  r <- terra::rasterize(terra::vect(xy[stats::complete.cases(xy), , drop = FALSE],
                                    type = "points", crs = ""),
                        r, fun = "count")
  r[is.na(r)] <- 0
  r <- r / sum(terra::values(r), na.rm = TRUE)

  ## circular padding: tile left and right copies for wrapping
  l_copy <- r
  r_copy <- r
  e <- terra::ext(r)
  terra::ext(l_copy) <- terra::ext(e[1] - 2 * pi, e[2] - 2 * pi, e[3], e[4])
  terra::ext(r_copy) <- terra::ext(e[1] + 2 * pi, e[2] + 2 * pi, e[3], e[4])
  r_padded <- terra::merge(l_copy, r, r_copy)

  ## bilinear interpolation to smooth
  target <- terra::rast(ncol = 150, nrow = 150,
                        xmin = -pi, xmax = pi,
                        ymin = 0, ymax = max_step,
                        crs = NA)
  r_smooth <- terra::resample(r_padded, target, method = "bilinear")
  r_smooth[r_smooth < 0] <- 0
  r_smooth <- r_smooth / sum(terra::values(r_smooth), na.rm = TRUE)

  r_smooth
}

#' Parametric step-turn probability using marginal distributions
#'
#' Fits a Weibull distribution to step lengths and a von Mises distribution
#' to turning angles, then returns the product of marginal densities
#' (independence copula).
#' @param step_length Numeric vector of step lengths.
#' @param turn_angle Numeric vector of turning angles (radians).
#' @return Numeric vector of probabilities (same length as input).
#' @noRd
.parametric_step_turn_prob <- function(step_length, turn_angle) {
  n <- length(step_length)
  prob <- rep(NA_real_, n)

  ## fit Weibull to positive step lengths
  valid_steps <- step_length[!is.na(step_length) & step_length > 0]
  if (length(valid_steps) < 3) return(prob)
  step_fit <- suppressWarnings(MASS::fitdistr(valid_steps, "weibull"))

  ## fit von Mises to turning angles
  valid_turns <- turn_angle[!is.na(turn_angle)]
  if (length(valid_turns) < 3) return(prob)
  circ_turns <- circular::circular(valid_turns, type = "angles",
                                   units = "radians")
  vm_fit <- circular::mle.vonmises(circ_turns)

  ## compute marginal densities
  step_dens <- rep(NA_real_, n)
  turn_dens <- rep(NA_real_, n)

  ok_s <- !is.na(step_length) & step_length > 0
  step_dens[ok_s] <- stats::dweibull(step_length[ok_s],
                                     shape = step_fit$estimate["shape"],
                                     scale = step_fit$estimate["scale"])

  ok_t <- !is.na(turn_angle)
  turn_dens[ok_t] <- as.numeric(circular::dvonmises(
    circular::circular(turn_angle[ok_t], type = "angles", units = "radians"),
    mu = vm_fit$mu, kappa = vm_fit$kappa))

  ## independence copula: product of marginals
  both_ok <- ok_s & ok_t
  prob[both_ok] <- step_dens[both_ok] * turn_dens[both_ok]
  prob
}

#' Optimise autodiff_alpha by maximising the trimmed mean of log-probabilities
#'
#' @param step_turn_prob Numeric vector of step-turn probabilities.
#' @param delta_step_prob Numeric vector of delta-step probabilities.
#' @param delta_turn_prob Numeric vector of delta-turn probabilities.
#' @return Numeric scalar: the optimal alpha.
#' @noRd
.optimize_alpha <- function(step_turn_prob, delta_step_prob, delta_turn_prob) {
  objective <- function(alpha, stp, dsp, dtp) {
    ad <- (dsp * dtp)^alpha
    joint <- stp * ad
    valid <- joint > 0 & !is.na(joint)
    if (sum(valid) < 3) return(Inf)
    -mean(log(joint[valid]), trim = 0.05)
  }
  result <- stats::optimize(objective, interval = c(0, 2),
                            stp = step_turn_prob,
                            dsp = delta_step_prob,
                            dtp = delta_turn_prob)
  result$minimum
}


#' Build prob-primitive reference distribution
#'
#' Constructs the (ref_step, ref_turn, ref_deltaStep, ref_deltaTurn,
#' ref_delta_gaps) list for a given reference move2 object.  Used by
#' both the legacy `reference =` path (one global ref) and the new
#' `pool_by =` path (one ref per pool group).
#'
#' @param reference move2 object to compute the reference from.
#' @param time_normalize Logical; if TRUE use speed / angular velocity.
#' @return A 5-element named list.
#' @noRd
.build_prob_ref_args <- function(reference, time_normalize) {
  ref_step <- as.numeric(move2::mt_distance(reference, units = "m"))
  ref_turn <- as.numeric(move2::mt_turnangle(reference))
  if (time_normalize) {
    ref_tl <- as.numeric(move2::mt_time_lags(reference, units = "secs"))
    ref_na_tl <- is.na(ref_tl)
    ref_step <- ifelse(ref_na_tl, NA_real_, ref_step / ref_tl)
    ref_turn <- ifelse(ref_na_tl, NA_real_, ref_turn / ref_tl)
    ref_delta_gaps <- c(NA, (ref_tl[-length(ref_tl)] + ref_tl[-1]) / 2)
  } else {
    ref_delta_gaps <- rep(1, length(ref_step))
  }
  list(ref_step       = ref_step,
       ref_turn       = ref_turn,
       ref_deltaStep  = c(NA, diff(ref_step)),
       ref_deltaTurn  = c(NA, diff(ref_turn)),
       ref_delta_gaps = ref_delta_gaps)
}


#' @param x A move2 object (original, full track).
#' @param iterations Maximum number of passes.
#' @param ... All other arguments forwarded to mt_flag_outliers (with iterations=1).
#' @return A move2 object with outlier flags mapped back to the original rows.
#' @noRd
.iterative_flag <- function(x, iterations, threshold, prob_type, remove, plot,
                            drop_na, autodiff_alpha, method,
                            quality_columns, time_normalize,
                            threshold_type, step_transform = "log",
                            step_floor = 0,
                            reference = NULL,
                            pool_by = NULL,
                            silent = FALSE) {
  say <- .say(silent)
  n <- nrow(x)
  all_outlier <- rep(FALSE, n)
  current <- x
  ## map from current row indices to original row indices
  current_to_orig <- seq_len(n)

  for (iter in seq_len(iterations)) {
    say(sprintf("--- Iteration %d of %d ---", iter, iterations))
    result <- mt_flag_outliers(current, threshold = threshold,
                               prob_type = prob_type, remove = FALSE,
                               plot = FALSE, drop_na = FALSE,
                               autodiff_alpha = autodiff_alpha,
                               method = method,
                               iterations = 1,
                               quality_columns = quality_columns,
                               time_normalize = time_normalize,
                               threshold_type = threshold_type,
                               step_transform = step_transform,
                               step_floor = step_floor,
                               reference = reference,
                               pool_by = pool_by,
                               silent = silent)
    new_outliers <- which(result$is_outlier)
    if (length(new_outliers) == 0) {
      say("No new outliers found. Stopping iteration.")
      break
    }
    ## map back to original indices
    orig_idx <- current_to_orig[new_outliers]
    all_outlier[orig_idx] <- TRUE
    say(sprintf("  Flagged %d new outliers (total: %d).",
                    length(new_outliers), sum(all_outlier)))

    ## mask outliers and update mapping
    keep <- !result$is_outlier
    current <- current[keep, ]
    current_to_orig <- current_to_orig[keep]

    if (nrow(current) < 4) {
      say("Too few locations remaining after masking. Stopping iteration.")
      break
    }
  }

  ## run one final pass on the original object (iterations=1) to get

  ## probabilities, then overwrite the is_outlier column with cumulative flags
  final <- mt_flag_outliers(x, threshold = threshold, prob_type = prob_type,
                            remove = FALSE, plot = FALSE, drop_na = FALSE,
                            autodiff_alpha = autodiff_alpha, method = method,
                            iterations = 1,
                            quality_columns = quality_columns,
                            time_normalize = time_normalize,
                            threshold_type = threshold_type,
                            step_transform = step_transform,
                            reference = reference,
                            pool_by = pool_by,
                            silent = silent)
  final$is_outlier <- all_outlier
  final <- .attach_prob_loglr(final)

  ## update outlier_percentile for the combined flags
  n_outliers <- sum(all_outlier)
  n_valid <- sum(!final$is_na_prob)
  prob_label <- switch(prob_type,
    joint = "Joint probability", step_turn = "Step-turn probability",
    delta_step = "Delta step probability",
    delta_turn = "Delta turn probability",
    custom = "Custom probability (step-turn x delta-step)")
  say(sprintf(
    "Iterative total: %d outliers (%.2f%% of %d locations) based on %s.",
    n_outliers, 100 * n_outliers / n_valid, n_valid, prob_label
  ))

  if (plot) {
    prob_col <- switch(prob_type,
      joint = "joint_prob", step_turn = "step_turn_prob",
      delta_step = "delta_step_prob", delta_turn = "delta_turn_prob",
      custom = "custom_prob")
    .plot_outliers(final, prob_col, prob_label, final$is_outlier,
                   final$is_na_prob, drop_na, n_outliers, sum(final$is_na_prob))
  }

  if (remove) {
    to_remove <- final$is_outlier
    if (drop_na) to_remove <- to_remove | final$is_na_prob
    return(final[!to_remove, ])
  }
  final
}


#' Two-panel diagnostic plot for outlier detection
#' @noRd
.plot_outliers <- function(x, prob_col, prob_label, is_outlier, is_na_prob,
                           drop_na, n_outliers, n_na) {
  message("Creating diagnostic plot...")

  coords <- sf::st_coordinates(x)
  log_prob <- log10(x[[prob_col]])

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))

  prob_colors <- grDevices::colorRampPalette(
    c("red", "yellow", "green", "blue")
  )(100)

  ## Panel 1: all locations coloured by probability
  graphics::plot(coords, type = "l", col = "grey80",
                 main = "Probability", xlab = "X", ylab = "Y", asp = 1)

  valid <- !is_na_prob
  if (any(valid)) {
    lp <- log_prob[valid]
    rng <- range(lp, na.rm = TRUE)
    if (diff(rng) > 0) {
      ci <- ceiling(((lp - rng[1]) / (rng[2] - rng[1])) * 99) + 1
    } else {
      ci <- rep(50, length(lp))
    }
    ## plot high-prob first, low-prob (outliers) on top
    ord <- order(lp, decreasing = TRUE)
    graphics::points(coords[valid, , drop = FALSE][ord, ],
                     col = prob_colors[ci[ord]], pch = 19, cex = 0.8)
  }
  if (any(is_na_prob)) {
    graphics::points(coords[is_na_prob, , drop = FALSE],
                     col = "grey50", pch = 19, cex = 0.8)
  }
  graphics::legend("topright",
                   legend = c("Low prob", "High prob", "NA"),
                   col = c(prob_colors[1], prob_colors[100], "grey50"),
                   pch = 19, cex = 0.8, bty = "n")

  ## Panel 2: kept vs removed
  to_remove <- is_outlier
  if (drop_na) to_remove <- to_remove | is_na_prob
  kept <- which(!to_remove)
  removed <- which(to_remove)

  graphics::plot(coords, type = "n", main = "Outliers removed",
                 xlab = "X", ylab = "Y", asp = 1)
  if (length(kept) > 1) {
    graphics::lines(coords[kept, ], col = "grey80")
  }
  if (length(kept) > 0) {
    graphics::points(coords[kept, , drop = FALSE],
                     col = "steelblue", pch = 19, cex = 0.6)
  }
  if (length(removed) > 0) {
    graphics::points(coords[removed, , drop = FALSE],
                     col = "red", pch = 19, cex = 1.0)
  }
  graphics::legend("topright",
                   legend = c("Kept", "Removed"),
                   col = c("steelblue", "red"),
                   pch = 19, cex = 0.8, bty = "n")

  ind_label <- paste(unique(move2::mt_track_id(x)), collapse = ", ")
  pct <- round(100 * length(removed) / nrow(x), 1)
  graphics::title(
    paste0(length(removed), " removed (", pct, "%) --", ind_label),
    outer = TRUE
  )
}


## Shared convention (DESIGN_evidence_accumulation.md): a detector's score
## -> a signed log-likelihood-ratio in nat units, zero-crossing at the
## detector's OWN flag boundary.  For the probability detector the score is
## the surprisal neglogp = -log(joint_prob) under the inlier model; the
## boundary is the least-surprising flagged fix.  Result > 0 exactly where
## the detector flags, < 0 below, NA where the detector abstains
## (invalid/NA probability).  Purely additive -- does not touch is_outlier.
## @keywords internal
## Shared: given a per-fix surprisal (neglogp, higher = more outlier, in
## nat units) and the detector's flag, return a signed log-LR that crosses
## zero at the detector's OWN flag boundary (the least-surprising flagged
## fix).  NA surprisal = abstain.
.loglr_from_surprisal <- function(neglogp, is_outlier) {
  neglogp[!is.finite(neglogp)] <- NA_real_
  flagged <- which((is_outlier %in% TRUE) & is.finite(neglogp))
  thr <- if (length(flagged)) {
    min(neglogp[flagged], na.rm = TRUE)
  } else {
    m <- suppressWarnings(max(neglogp, na.rm = TRUE))
    if (is.finite(m)) m + 1e-9 else 0
  }
  neglogp - thr
}

## Group-aware loglr: apply .loglr_from_surprisal within each group so the
## flag boundary is computed per track / per pool, never across them.
## Uses dplyr group_by (preserves row order under mutate).  group = NULL
## reduces to the global (single-track) computation.
## @keywords internal
.loglr_grouped <- function(neglogp, is_outlier, group = NULL) {
  if (is.null(group) || length(unique(group)) <= 1L) {
    return(.loglr_from_surprisal(neglogp, is_outlier))
  }
  d <- data.frame(np = neglogp, io = is_outlier %in% TRUE,
                  .grp = as.character(group), stringsAsFactors = FALSE)
  d <- dplyr::group_by(d, dplyr::across(dplyr::all_of(".grp")))
  d <- dplyr::mutate(d, ll = .loglr_from_surprisal(.data[["np"]],
                                                   .data[["io"]]))
  dplyr::pull(dplyr::ungroup(d), "ll")
}

## Probability detector: surprisal = -log(joint_prob).
.attach_prob_loglr <- function(obj) {
  ## flag-source column for votability (== is_outlier for a single detector).
  obj$flagged_by_prob <- obj$is_outlier
  if (is.null(obj$joint_prob)) return(obj)
  neglogp <- -log(obj$joint_prob)
  if (!is.null(obj$is_na_prob)) {
    neglogp[as.logical(obj$is_na_prob) %in% TRUE] <- NA_real_
  }
  obj$loglr_prob <- .loglr_grouped(neglogp, obj$is_outlier,
                                   move2::mt_track_id(obj))
  obj
}

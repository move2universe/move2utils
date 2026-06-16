## Shared gap-threshold helper for outlier detection.
##
## Extracted from mt_flag_outliers() so the bridge-residual detector
## can apply the same broken-stick + tail-decay inflection logic to
## a different underlying score.  Kept internal (dot-prefixed).


#' Broken-stick plus tail-decay inflection threshold on a numeric vector
#'
#' Given a vector of log-scores where \emph{low values indicate outliers}
#' (e.g. log joint probability, or \code{-log(bridge_eta)}), identify
#' a natural break separating an outlier tail from the bulk.
#'
#' @param log_scores Numeric vector, typically log-probabilities or
#'   \code{-log(score)}.  NA values are tolerated and excluded.
#' @param threshold Multiplier controlling the strictness of the break
#'   detection.  A gap must exceed \code{threshold} times the local
#'   expectation (broken-stick) or a noise floor (inflection) to count.
#'   Default 3 is HEURISTIC -- the "3-sigma rule" applied loosely to
#'   the broken-stick gap distribution. Plausible range 2--5; lower
#'   values are more sensitive (more candidate breaks accepted),
#'   higher values are more conservative.
#' @param search_frac Fraction of the sorted values from the lower tail
#'   to consider as candidate break locations.  Default \code{1/7} is
#'   HEURISTIC (~14\% of the sorted values, the conventional "search the
#'   lowest sextile to find the outlier-tail break"). Plausible range:
#'   1/10 to 1/5.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{is_outlier}}{Logical vector of the same length as
#'       \code{log_scores}, \code{TRUE} for values below the break.}
#'     \item{\code{break_value}}{Numeric.  Values strictly less than
#'       this are outliers.  \code{NA} if no break found.}
#'     \item{\code{method}}{Character: \code{"inflection"},
#'       \code{"broken_stick"}, \code{"none"}, or \code{"insufficient"}.}
#'     \item{\code{percentile}}{Numeric vector of empirical percentiles
#'       (0--100) for all input values, for diagnostic purposes.}
#'   }
#'
#' @keywords internal
.gap_threshold_lower <- function(log_scores, threshold = 3,
                                  search_frac = 1 / 7) {
  n <- length(log_scores)
  is_na  <- is.na(log_scores)
  is_fin <- !is_na & is.finite(log_scores)
  valid  <- log_scores[is_fin]

  ## percentile diagnostic (always returned)
  pct <- rep(NA_real_, n)
  if (length(valid) >= 2) {
    pct[is_fin] <- round(100 * stats::ecdf(valid)(valid), 4)
  }

  ## insufficient data
  if (length(valid) < 10) {
    return(list(is_outlier = rep(FALSE, n),
                break_value = NA_real_,
                method = "insufficient",
                percentile = pct))
  }

  sorted   <- sort(valid)
  n_sorted <- length(sorted)
  gaps     <- diff(sorted)
  n_gaps   <- length(gaps)

  ## not enough distinct values
  if (n_gaps < 2 || (sorted[n_sorted] - sorted[1]) == 0) {
    return(list(is_outlier = rep(FALSE, n),
                break_value = NA_real_,
                method = "none",
                percentile = pct))
  }

  ## --- stage 1: broken-stick screening ---
  total_range  <- sorted[n_sorted] - sorted[1]
  harmonic     <- rev(cumsum(1 / seq(n_gaps, 1)))
  expected_gaps <- (total_range / n_gaps) * harmonic / sum(harmonic) * n_gaps
  gap_ratio    <- gaps / pmax(expected_gaps, .Machine$double.eps)

  search_n      <- max(2, n_sorted %/% round(1 / search_frac))
  search_region <- seq_len(min(search_n, n_gaps))

  ## --- stage 2: tail-decay inflection ---
  bulk_start  <- sorted[min(search_n + 1, n_sorted)]
  tail_length <- bulk_start - sorted[seq_len(search_n)]

  inflection_idx <- 0L
  if (length(tail_length) >= 4) {
    d2 <- diff(diff(tail_length))
    right_half <- ceiling(length(d2) / 2):length(d2)
    noise_floor <- if (length(right_half) > 0)
                     stats::median(abs(d2[right_half]))
                   else stats::median(abs(d2))
    if (is.finite(noise_floor) && noise_floor > 0) {
      below_noise <- which(abs(d2) < threshold * noise_floor)
      if (length(below_noise) > 0) {
        inflection_idx <- below_noise[1] + 1L
      }
    }
  }

  ## --- combine ---
  is_outlier  <- rep(FALSE, n)
  break_value <- NA_real_
  method      <- "none"

  if (inflection_idx > 1L && inflection_idx <= search_n) {
    break_value <- sorted[inflection_idx]
    method      <- "inflection"
  } else {
    ratios_in_region <- gap_ratio[search_region]
    best_idx  <- search_region[which.max(ratios_in_region)]
    best_ratio <- gap_ratio[best_idx]
    if (is.finite(best_ratio) && best_ratio > threshold) {
      break_value <- sorted[best_idx]
      method      <- "broken_stick"
    }
  }

  if (!is.na(break_value)) {
    is_outlier[is_fin] <- log_scores[is_fin] < break_value
  }

  list(is_outlier = is_outlier,
       break_value = break_value,
       method = method,
       percentile = pct)
}


#' Entropy-valley threshold on a numeric vector
#'
#' Given a vector of log-scores where \emph{low values indicate outliers},
#' estimate the density via KDE and search for the deepest local minimum
#' (valley) below the main mode.  A valley indicates a natural separation
#' between an outlier regime and the bulk.  If no valley meets the depth
#' criterion, no outliers are declared -- so this threshold is safe on
#' clean data where the gap detector over-flags.
#'
#' @param log_scores Numeric vector.  NAs are tolerated.
#' @param threshold Numeric in (0, 1).  Maximum allowed ratio of valley
#'   density to peak density.  A valley whose density exceeds
#'   \code{threshold * peak_density} is not deep enough to justify
#'   splitting the distribution.  Default 0.3 is the unified package-
#'   wide entropy default, validated by the 2026-05-06 Raven sensitivity
#'   sweep across 65 stratified Movebank tracks (the only level inside
#'   the strict stability window for cohort flag rate; K-W p = 0.076,
#'   i.e. cohort flag rate is statistically insensitive in the
#'   0.3--0.7 range).  Plausible range: 0.3--0.7.
#' @param n_grid KDE grid size, default 512 (standard density-estimation
#'   grid; grid resolution rarely matters for valley detection).
#'
#' @return A list with \code{is_outlier}, \code{break_value},
#'   \code{method} (one of \code{"valley"}, \code{"none"},
#'   \code{"insufficient"}), and \code{percentile}.
#'
#' @keywords internal
.entropy_threshold_lower <- function(log_scores, threshold = 0.3,
                                       n_grid = 512) {
  n <- length(log_scores)
  is_na  <- is.na(log_scores)
  is_fin <- !is_na & is.finite(log_scores)
  valid  <- log_scores[is_fin]

  pct <- rep(NA_real_, n)
  if (length(valid) >= 2) {
    pct[is_fin] <- round(100 * stats::ecdf(valid)(valid), 4)
  }

  if (length(valid) < 10) {
    return(list(is_outlier = rep(FALSE, n),
                break_value = NA_real_,
                method = "insufficient",
                percentile = pct))
  }

  d <- stats::density(valid, n = n_grid)
  peak_idx <- which.max(d$y)
  peak_density <- d$y[peak_idx]
  peak_x <- d$x[peak_idx]

  ## search for valleys below (left of) the peak
  lower <- which(d$x < peak_x)
  valley_x <- NA_real_
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

  is_outlier <- rep(FALSE, n)
  break_value <- NA_real_
  method <- "none"

  if (!is.na(valley_x) && valley_depth < threshold) {
    break_value <- valley_x
    method <- "valley"
    is_outlier[is_fin] <- log_scores[is_fin] < valley_x
  }

  list(is_outlier = is_outlier,
       break_value = break_value,
       method = method,
       percentile = pct)
}


#' Dip-test-validated combined threshold
#'
#' Runs the entropy detector first (strict; returns no flags on a
#' clean unimodal distribution).  If entropy finds no valley, falls
#' back to the gap/broken-stick detector, but accepts its proposed
#' break only if Hartigan's dip test indicates the distribution is
#' significantly multimodal (\code{p < dip_alpha}).  This couples
#' gap's sensitivity with a formal statistical check on the
#' unimodality null, so gap-over-flagging on legitimate heavy tails
#' is suppressed.
#'
#' @param log_scores Numeric vector where low values indicate outliers.
#' @param entropy_threshold Density-ratio for the entropy detector.
#'   \code{NULL} (default) defers to
#'   \code{\link{.entropy_threshold_lower}}'s own formal default --
#'   the single source of truth for the package-wide entropy default.
#'   See that function for the heuristic basis and plausible range.
#' @param gap_threshold Break-size multiplier for the gap detector.
#'   \code{NULL} (default) defers to
#'   \code{\link{.gap_threshold_lower}}'s own formal default.
#' @param dip_alpha Significance threshold for the dip test on the
#'   unimodality null.  Default 0.05 is the standard significance
#'   level (Fisher convention); not derived from the package's
#'   benchmark data.
#'
#' @return Same structure as \code{\link{.gap_threshold_lower}} with
#'   an added \code{dip_p} entry when the gap fallback was invoked.
#'
#' @keywords internal
.entropy_or_dip_gap_threshold_lower <- function(
    log_scores,
    entropy_threshold = NULL,
    gap_threshold     = NULL,
    dip_alpha         = 0.05) {

  ## NULL defers to the leaf's own formal default (single source of
  ## truth); a non-NULL value is forwarded explicitly.  Pattern A from
  ## audits/2026-05-25-parameter-propagation/findings.md §1.
  e_args <- list(log_scores)
  if (!is.null(entropy_threshold)) e_args$threshold <- entropy_threshold
  e <- do.call(.entropy_threshold_lower, e_args)
  if (!is.na(e$break_value)) return(c(e, list(dip_p = NA_real_)))

  g_args <- list(log_scores)
  if (!is.null(gap_threshold)) g_args$threshold <- gap_threshold
  g <- do.call(.gap_threshold_lower, g_args)
  if (is.na(g$break_value)) return(c(g, list(dip_p = NA_real_)))

  ## Gap proposed a candidate. Validate against unimodality null
  ## via Hartigan's dip test. Only accept if the distribution is
  ## significantly multimodal -- otherwise the candidate is almost
  ## certainly a heavy-tail artifact.
  valid <- log_scores[is.finite(log_scores)]
  if (length(valid) < 10 || !requireNamespace("diptest", quietly = TRUE)) {
    ## Without diptest we cannot validate; reject gap to stay conservative.
    g$is_outlier[]   <- FALSE
    g$break_value    <- NA_real_
    g$method         <- "gap_rejected_no_diptest"
    return(c(g, list(dip_p = NA_real_)))
  }

  dip <- diptest::dip.test(valid)
  if (dip$p.value < dip_alpha) {
    g$method <- paste0(g$method, "_dip_validated")
    return(c(g, list(dip_p = dip$p.value)))
  }

  ## Not significantly multimodal -- reject gap's candidate.
  g$is_outlier[]   <- FALSE
  g$break_value    <- NA_real_
  g$method         <- "gap_rejected_by_dip"
  c(g, list(dip_p = dip$p.value))
}

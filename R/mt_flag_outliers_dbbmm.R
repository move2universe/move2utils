## Isotropic special case of the state-aware bridge primitive.
##
## mt_flag_outliers_dbgb() is the canonical state-aware bridge
## primitive: it estimates per-axis (parallel / orthogonal) motion
## variance via dBGB and runs a per-channel envelope rule on the
## resulting Z scores.  Under the constraint sigma_para = sigma_orth
## the dBGB model reduces to dBBMM, the residual decomposition
## collapses to its Euclidean magnitude, and the envelope's per-axis
## tests become inert -- the chisq channel is the only one that fires.
## That submodel is the historical "dBBMM" primitive.
##
## This wrapper preserves the historical entry point and gives users
## a one-knob "go faster, give up the directional decomposition"
## option.  It delegates to mt_flag_outliers_dbgb() with variance
## fixed at "dbbmm".


#' Flag outliers via dBBMM-Z (isotropic submodel of dBGB)
#'
#' Convenience wrapper around \code{\link{mt_flag_outliers_dbgb}} that
#' fixes \code{variance = "dbbmm"}.  Mathematically the
#' \eqn{\sigma^2_{\parallel} = \sigma^2_{\perp}} submodel of the
#' dBGB primitive: the bivariate residual is normalised by a single
#' isotropic motion variance and the per-axis decomposition is
#' discarded.  Faster (one 1-D Brent fit per window rather than two)
#' and emits a simpler output schema (\code{sigma2_motion} instead of
#' \code{sigma2_motion_para} / \code{sigma2_motion_orth};
#' \code{bridge_z_para} / \code{bridge_z_orth} not present).  The
#' flagging Z is the chisq scalar Rayleigh-Z and the parametric
#' threshold helpers (Bonferroni / BH-FDR) apply directly.
#'
#' @param x A \code{move2} object.  Single- or multi-track.  Auto-
#'   projected to a per-track local AEQD if input is lon/lat.
#' @param ... Further arguments passed to
#'   \code{\link{mt_flag_outliers_dbgb}}.  Do not pass \code{variance};
#'   it is fixed at \code{"dbbmm"} here -- call
#'   \code{\link{mt_flag_outliers_dbgb}} directly to vary it.  The
#'   per-channel \code{z_threshold_method} options
#'   (\code{"bonferroni"}, \code{"bh_fdr"}, \code{"gap"}) are
#'   forwarded; under \code{variance = "dbbmm"} only the chisq
#'   channel exists, so the envelope rule reduces to a scalar Z test.
#'
#' @return The input object with bridge-Z outlier-detection columns
#'   attached; see \code{\link{mt_flag_outliers_dbgb}} for the column
#'   schema.  Per-axis columns
#'   (\code{bridge_z_para}, \code{bridge_z_orth},
#'   \code{sigma2_motion_para}, \code{sigma2_motion_orth}) are not
#'   emitted under \code{variance = "dbbmm"}.
#'
#' @seealso \code{\link{mt_flag_outliers_dbgb}} for the full state-
#'   aware bridge primitive with directional decomposition;
#'   \code{\link{mt_flag_outliers_bridge}} for the leverage-immune
#'   static bridge primitive.
#'
#' @examples
#' \dontrun{
#' library(move2)
#' x <- mt_read(system.file("extdata/Pettstadt1-14053.csv.gz",
#'                            package = "move2utils"))
#' x <- x[!sf::st_is_empty(x), ]
#' res <- mt_flag_outliers_dbbmm(x)
#' }
#'
#' @export
mt_flag_outliers_dbbmm <- function(x, ...) {
  dots <- list(...)
  if ("variance" %in% names(dots))
    rlang::abort(paste0(
      "`variance` cannot be set in mt_flag_outliers_dbbmm; ",
      "it is fixed at \"dbbmm\".  Call mt_flag_outliers_dbgb() ",
      "directly to vary it."),
      class = "move2utils_mt_flag_outliers_dbbmm_fixed_variance")
  mt_flag_outliers_dbgb(x, variance = "dbbmm", ...)
}

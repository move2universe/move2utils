#' Allometric estimate of physiological maximum speed
#'
#' Returns a body-mass / locomotor-mode prediction of the maximum
#' burst speed an animal of the given mass and mode is physically
#' capable of, using the general scaling law of Hirt et al. (2017).
#' Intended as a principled default for the \code{v_max} cap used by
#' \code{\link{mt_clean_track}} and \code{\link{mt_peel_speed}} when
#' the user does not have a species-specific number.
#'
#' @details
#' The Hirt et al. (2017) model is a time-dependent saturation of a
#' power-law scaling of theoretical maximum speed with body mass:
#' \deqn{v_{max} = a M^{b} \, (1 - e^{-h M^{i}})}
#' where the first factor is the theoretical aerobic ceiling and the
#' saturation term captures the finite anaerobic energy budget that
#' limits large animals to a fraction of that ceiling.  The fitted
#' parameters per locomotion mode are taken directly from
#' Supplementary Table 4 of Hirt et al. (2017) (M in kg, output v in
#' km/h, internally converted to m/s):
#' \tabular{lllll}{
#'   mode      \tab a            \tab b           \tab h           \tab i \cr
#'   flying    \tab 142.8 +- 16.7\tab 0.24 +- 0.01 \tab 2.4 +- 1.4 \tab -0.72 +- 0.26 \cr
#'   running   \tab  25.5 +- 0.84\tab 0.26 +- 0.006\tab  22 +- 7.6 \tab -0.6  +- 0.05 \cr
#'   swimming  \tab  11.2 +- 0.91\tab 0.36 +- 0.02 \tab 19.5 +- 13.6\tab -0.56 +- 0.07
#' }
#' Reported predictive accuracy R^2 = 0.893 across 622 data points
#' from 474 species spanning 3e-8 to 1.084e5 kg.
#'
#' \strong{Scope and caveats.}
#' \itemize{
#'   \item Hirt et al. excluded vertical (gravity-assisted) speeds
#'         from their dataset.  This is the right scope for
#'         \code{move2utils}, which operates on horizontally projected
#'         GPS / Argos / satellite tracking data.  A peregrine stoop
#'         at 89 m/s in 3D projects to a much smaller horizontal
#'         step speed and falls under the model's prediction.  For
#'         3D-tracked diving / stooping data, the user should
#'         override with species-specific aerodynamic literature.
#'   \item The model's data are predominantly maximum anaerobic burst
#'         speeds -- the impossibility ceiling, not the cruise speed.
#'         This is the right semantic for \code{v_max}: outliers are
#'         transitions an animal cannot physically perform.
#'   \item Locomotor specialists (cheetah ~29 m/s, pronghorn antelope
#'         ~26 m/s, sailfish in water) sit in the upper tail of the
#'         residual distribution and can exceed the central
#'         prediction.  When a published species-specific maximum is
#'         available, use it directly instead of the allometric
#'         default.
#'   \item The returned prediction interval is propagated from the
#'         fitted parameter standard errors via the delta method.  It
#'         describes parameter uncertainty in the model fit, not the
#'         species-to-species residual scatter; the latter implies a
#'         further roughly factor-of-2 spread documented in the
#'         original paper as the model R^2 = 0.893.
#' }
#'
#' @param mass Numeric scalar.  Body mass in kilograms (kg).  Must be
#'   positive.  A warning is issued when \code{mass} lies outside the
#'   range of the Hirt et al. (2017) dataset (3e-8 to 1.084e5 kg).
#' @param mode Character.  One of \code{"flying"}, \code{"running"},
#'   \code{"swimming"}.
#' @param ci_level Numeric in (0, 1).  Confidence level for the
#'   parameter-uncertainty interval reported as the \code{"ci"}
#'   attribute.  Default 0.95.
#'
#' @return A length-1 numeric vector containing the central prediction
#'   in m/s.  Attached attributes:
#'   \describe{
#'     \item{\code{ci}}{Length-2 numeric, lower and upper bound of the
#'       parameter-uncertainty interval at \code{ci_level}, in m/s.}
#'     \item{\code{kmh}}{Length-1 numeric, central prediction in km/h
#'       (the original Hirt unit).}
#'     \item{\code{mass}}{The supplied body mass.}
#'     \item{\code{mode}}{The supplied locomotor mode.}
#'     \item{\code{reference}}{A string citing Hirt et al. (2017).}
#'   }
#'   Pass directly to \code{\link{mt_clean_track}} via the
#'   \code{v_max} argument; numeric coercion strips the attributes
#'   and yields the central estimate.  \code{mt_clean_track()} also
#'   accepts the \code{(mass, mode)} pair directly and runs the
#'   estimator internally; that is the recommended path.
#'
#' @examples
#' ## golden eagle, ~5 kg
#' v_phys_estimate(5, "flying")
#'
#' ## red fox, ~6 kg
#' v_phys_estimate(6, "running")
#'
#' ## bottlenose dolphin, ~250 kg
#' v_phys_estimate(250, "swimming")
#'
#' \dontrun{
#' ## use as the principled v_max default
#' clean <- mt_clean_track(track,
#'                          v_max = v_phys_estimate(mass = 5,
#'                                                      mode = "flying"))
#'
#' ## equivalent and more idiomatic: pass (mass, mode) directly
#' clean <- mt_clean_track(track, mass = 5, mode = "flying")
#' }
#'
#' @references
#' Hirt, M. R., Jetz, W., Rall, B. C., Brose, U. (2017).  A general
#' scaling law reveals why the largest animals are not the fastest.
#' \emph{Nature Ecology & Evolution} 1, 1116-1122.
#' \doi{10.1038/s41559-017-0241-4}
#'
#' @seealso \code{\link{mt_clean_track}}, \code{\link{mt_peel_speed}},
#'   \code{\link{mt_suggest_speed_cap}}
#' @export
v_phys_estimate <- function(mass, mode = c("flying", "running", "swimming"),
                                ci_level = 0.95) {

  ## ---- input validation ------------------------------------------
  mode <- match.arg(mode)
  if (!is.numeric(mass) || length(mass) != 1L || is.na(mass) || mass <= 0) {
    rlang::abort("`mass` must be a positive scalar (kg).",
                 class = "move2utils_v_phys_estimate_bad_mass")
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      ci_level <= 0 || ci_level >= 1) {
    rlang::abort("`ci_level` must be a scalar in (0, 1).",
                 class = "move2utils_v_phys_estimate_bad_ci_level")
  }
  hirt_range <- c(3e-8, 1.084e5)
  if (mass < hirt_range[1] || mass > hirt_range[2]) {
    rlang::warn(
      sprintf("mass = %g kg is outside the range of the Hirt et al. (2017) dataset (%g, %g) kg; the prediction is an extrapolation.",
              mass, hirt_range[1], hirt_range[2]),
      class = "move2utils_v_phys_estimate_extrapolation")
  }

  ## ---- Hirt et al. 2017 Supplementary Table 4 --------------------
  par_table <- list(
    flying   = list(a = 142.8, a_se = 16.7,
                    b = 0.24,  b_se = 0.01,
                    h = 2.4,   h_se = 1.4,
                    i = -0.72, i_se = 0.26),
    running  = list(a = 25.5,  a_se = 0.84,
                    b = 0.26,  b_se = 0.006,
                    h = 22,    h_se = 7.6,
                    i = -0.6,  i_se = 0.05),
    swimming = list(a = 11.2,  a_se = 0.91,
                    b = 0.36,  b_se = 0.02,
                    h = 19.5,  h_se = 13.6,
                    i = -0.56, i_se = 0.07)
  )
  p <- par_table[[mode]]

  ## ---- central prediction ----------------------------------------
  M  <- mass
  v_kmh <- p$a * M^p$b * (1 - exp(-p$h * M^p$i))
  v_ms  <- v_kmh / 3.6

  ## ---- delta-method parameter CI in km/h -------------------------
  ## v = a * M^b * (1 - exp(-h * M^i))
  ## let q = M^b, r = M^i, s = exp(-h * r)
  ## dv/da = q * (1 - s)
  ## dv/db = a * q * log(M) * (1 - s)
  ## dv/dh = a * q * s * r
  ## dv/di = a * q * s * h * r * log(M)
  q <- M^p$b
  r <- M^p$i
  s <- exp(-p$h * r)
  dv_da <- q * (1 - s)
  dv_db <- p$a * q * log(M) * (1 - s)
  dv_dh <- p$a * q * s * r
  dv_di <- p$a * q * s * p$h * r * log(M)

  var_v_kmh <- (dv_da * p$a_se)^2 + (dv_db * p$b_se)^2 +
               (dv_dh * p$h_se)^2 + (dv_di * p$i_se)^2
  se_v_kmh  <- sqrt(var_v_kmh)
  z         <- stats::qnorm(0.5 + ci_level / 2)
  ci_kmh    <- v_kmh + c(-1, 1) * z * se_v_kmh
  ci_kmh    <- pmax(ci_kmh, 0)
  ci_ms     <- ci_kmh / 3.6

  ## ---- assemble return value -------------------------------------
  out <- v_ms
  attr(out, "ci")        <- ci_ms
  attr(out, "kmh")       <- v_kmh
  attr(out, "mass")      <- mass
  attr(out, "mode")      <- mode
  attr(out, "reference") <- "Hirt et al. (2017) Nat. Ecol. Evol., doi:10.1038/s41559-017-0241-4"
  attr(out, "ci_level")  <- ci_level
  class(out) <- c("v_phys_estimate", class(out))
  out
}


#' @export
print.v_phys_estimate <- function(x, ...) {
  ci <- attr(x, "ci")
  cat(sprintf(
    "v_max (allometric, Hirt et al. 2017): %.2f m/s  (%.2f km/h)\n",
    as.numeric(x), attr(x, "kmh")))
  cat(sprintf(
    "  %.0f%% parameter CI: [%.2f, %.2f] m/s\n",
    100 * attr(x, "ci_level"), ci[1], ci[2]))
  cat(sprintf("  inputs:  mass = %g kg, mode = %s\n",
              attr(x, "mass"), attr(x, "mode")))
  cat(sprintf("  source:  %s\n", attr(x, "reference")))
  invisible(x)
}

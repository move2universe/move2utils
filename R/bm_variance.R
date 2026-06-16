#' Brownian Motion Variance Estimation
#'
#' Estimate the Brownian motion variance parameter from a set of locations
#' using a leave-one-out likelihood approach.
#'
#' @param time_lag Numeric vector of time lags (in minutes).
#' @param location_error Numeric vector of location errors per position.
#' @param x Numeric vector of x coordinates.
#' @param y Numeric vector of y coordinates.
#'
#' @return A list with `BMvar` (the estimated variance) and `cll` (the
#'   conditional log-likelihood).
#'
#' @details This is the pure R reference implementation, superseded by
#'   the C version in `bm_variance_c.c` for production use. It implements
#'   the leave-one-out approach of Horne et al. (2007).
#'
#' @references
#' Horne, J. S., Garton, E. O., Krone, S. M., & Lewis, J. S. (2007).
#' Analyzing animal movements using Brownian bridges. *Ecology*, 88(9),
#' 2354-2363.
#'
#' @keywords internal
.bm_variance <- function(time_lag, location_error, x, y) {
  n_locs <- unique(c(length(time_lag), length(location_error), length(x), length(y)))
  if (length(n_locs) != 1) {
    rlang::abort("Not an equal number of locations in call to .bm_variance",
                 class = "move2utils_internal_bm_variance_unequal_lengths")
  }

  T_jump <- alpha <- ztz <- loc_err_1 <- loc_err_2 <- NULL
  i <- 2
  while (i < n_locs) {
    t <- time_lag[i - 1] + time_lag[i]
    T_jump <- c(T_jump, t)
    a <- time_lag[i - 1] / t
    alpha <- c(alpha, a)
    u <- c(x[i - 1], y[i - 1]) + a * (c(x[i + 1], y[i + 1]) - c(x[i - 1], y[i - 1]))
    ztz <- c(ztz, (c(x[i], y[i]) - u) %*% (c(x[i], y[i]) - u))
    loc_err_1 <- c(loc_err_1, location_error[i - 1])
    loc_err_2 <- c(loc_err_2, location_error[i + 1])
    i <- i + 2
  }

  likelihood <- function(var, T_jump, alpha, loc_err_1, loc_err_2, ztz) {
    v <- T_jump * alpha * (1 - alpha) * var +
      ((1 - alpha)^2) * (loc_err_1^2) +
      (alpha^2) * (loc_err_2^2)
    l <- log((1 / (2 * pi * v)) * exp(-ztz / (2 * v)))
    if (any(is.na(l))) {
      rlang::abort("Internal error in variance estimation. Contact maintainer.",
                   class = "move2utils_internal_bm_variance_likelihood_na")
    }
    return(-sum(l))
  }

  lower <- 0
  upper <- 1e+15
  BMvar <- optimize(likelihood,
    lower = lower, upper = upper,
    T_jump = T_jump, alpha = alpha,
    loc_err_1 = loc_err_1, loc_err_2 = loc_err_2,
    ztz = ztz
  )

  if (any(BMvar$minimum %in% c(lower, upper))) {
    rlang::abort("Optimization failed! Consider changing map units.",
                 class = "move2utils_bm_variance_optim_failed")
  }

  if ((length(x) %% 2) == 0) {
    rlang::warn(
      "Even number of locations in variance function; last location unused.",
      class = "move2utils_bm_variance_even_n")
  }

  list(BMvar = BMvar$minimum, cll = -BMvar$objective)
}

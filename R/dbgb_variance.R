#' Parallel and Orthogonal Decomposition
#'
#' Decompose the displacement between a point and a reference position
#' into parallel and orthogonal components relative to a direction vector.
#'
#' @param mu Matrix of expected positions (n x 2).
#' @param direction_point Matrix or numeric vector defining the direction.
#' @param point Matrix of actual positions (n x 2).
#'
#' @return A matrix with columns `deltaPara` and `deltaOrth`.
#' @keywords internal
.delta_para_orth <- function(mu, direction_point, point) {
  if (!is.matrix(mu)) mu <- matrix(mu, ncol = 2, nrow = nrow(point), byrow = TRUE)
  if (!is.matrix(direction_point)) {
    direction_point <- matrix(direction_point, ncol = 2, nrow = nrow(point), byrow = TRUE)
  }

  C <- sqrt(rowSums((mu - direction_point)^2))
  B <- sqrt(rowSums((point - direction_point)^2))
  A <- sqrt(rowSums((point - mu)^2))

  tmp <- ((A^2 + C^2 - B^2) / (2 * A * C))

  same_loc <- (rowSums(mu == point) == ncol(mu))
  same_loc_dir <- (rowSums(mu == direction_point) == ncol(mu))

  tmp[same_loc_dir] <- sqrt(0.5)
  if (any(same_loc_dir)) {
    rlang::warn(
      "Brownian motion assumed, because no direction could be calculated",
      class = "move2utils_dbgb_variance_no_direction")
  }
  tmp[same_loc] <- 0
  tmp[tmp > 1] <- 1
  tmp[tmp < -1] <- -1

  theta <- acos(tmp)
  cbind(deltaPara = A * cos(theta), deltaOrth = A * sin(theta))
}


#' BGB Variance with Breakpoint Detection (C-accelerated)
#'
#' Per-window dBGB variance + breakpoint search.  Sequential search:
#' (1) no break, (2) best para-only break, (3) best orth-only break,
#' (4) cross-conditional refinement.  The whole search runs in C
#' (`bgb_var_window_c` in `src/bgb_var_window_c.c`); this wrapper
#' just prepares the inputs and re-shapes the output as the
#' historical R reference did.
#'
#' @details
#' The earlier all-R implementation drove `optim()` once per
#' breakpoint candidate and per axis combination, paying ~30-40 R-
#' level optim() calls per window with hundreds of `llBGBvar`
#' R<->C round-trips inside each one.  On a 1748-fix track that
#' scaled to ~200x slower than the dBBMM C kernel.  The C kernel
#' precomputes the leave-one-out quantities once per window, then
#' runs a 1D Brent's-method optimisation per axis (axes are
#' independent given fixed breaks under the diagonal-Sigma dBGB
#' model) and tests every breakpoint candidate inside the same
#' .Call.  Numerical agreement with the historical R reference is
#' verified by the test
#' \code{tests/testthat/test-bgb_var_break_c_vs_r.R}.
#'
#' The historical R reference is preserved as
#' \code{.bgb_var_break_r_reference()} for that numerical-
#' equivalence test; it is not used in production.
#'
#' @keywords internal
.bgb_var_break <- function(x_coords, y_coords, time_mins, location_error, margin) {
  n <- length(x_coords)
  if (length(location_error) == 1L) location_error <- rep(location_error, n)

  potential_breaks <- 2:(n - 1)
  margin_breaks <- potential_breaks[potential_breaks >= margin &
                                      potential_breaks <= (1 + n - margin) &
                                      (potential_breaks %% 2) == 1]
  if (length(margin_breaks) == 0L) {
    return(data.frame(paraSd = rep(NA_real_, n),
                      orthSd = rep(NA_real_, n)))
  }

  ## time_lag[i] = time_mins[i + 1] - time_mins[i]; the C kernel uses
  ## the same lag convention as the dBBMM C kernel and the R
  ## reference implementation, both of which sum consecutive lags to
  ## form the bridge T_jump.
  time_lag <- c(diff(time_mins), 0)

  res <- .Call("bgb_var_window_c",
               as.double(x_coords), as.double(y_coords),
               as.double(time_lag), as.double(location_error),
               as.integer(margin_breaks))

  data.frame(paraSd = res$paraSd, orthSd = res$orthSd)
}

#' Historical pure-R reference for .bgb_var_break (kept for numerical
#' equivalence testing)
#'
#' This is the implementation that shipped before the C kernel
#' \code{bgb_var_window_c} was introduced.  It uses optim()
#' (L-BFGS-B) once per candidate breakpoint and is ~200x slower than
#' the C path on a 1748-fix track.  Production code calls
#' \code{.bgb_var_break} (the C wrapper); this function exists solely
#' so that the equivalence test can compare per-axis sigmas against
#' an independent reference.
#'
#' @keywords internal
.bgb_var_break_r_reference <- function(x_coords, y_coords, time_mins,
                                       location_error, margin) {
  n <- length(x_coords)
  if (length(location_error) == 1L) location_error <- rep(location_error, n)
  coords <- cbind(x_coords, y_coords)

  is <- (1:n)[(1:n) %% 2 == 0]
  alphas <- (time_mins[is] - time_mins[is - 1]) /
    (time_mins[is + 1] - time_mins[is - 1])
  mus <- coords[is - 1, , drop = FALSE] +
    alphas * (coords[is + 1, , drop = FALSE] - coords[is - 1, , drop = FALSE])
  para_orth <- .delta_para_orth(mus, coords[is + 1, , drop = FALSE],
                                coords[is, , drop = FALSE])

  errs <- alphas^2 * location_error[is + 1]^2 + (1 - alphas)^2 * location_error[is - 1]^2
  sd_mul <- alphas * (1 - alphas) * (time_mins[is + 1] - time_mins[is - 1])
  n_pairs <- length(errs)

  potential_breaks <- 2:(n - 1)
  margin_breaks <- potential_breaks[potential_breaks >= margin &
                                      potential_breaks <= (1 + n - margin) &
                                      (potential_breaks %% 2) == 1]

  eval_ll <- function(para_sd_vec, orth_sd_vec) {
    sp <- sqrt(errs + para_sd_vec^2 * sd_mul)
    so <- sqrt(errs + orth_sd_vec^2 * sd_mul)
    .Call("llBGBvar", cbind(sp, so)^2, para_orth)
  }

  optimize_sds <- function(para_break, orth_break) {
    init <- c(paraBefore = 100, orthBefore = 100)
    if (!is.na(para_break)) init <- c(init, paraAfter = 100)
    if (!is.na(orth_break)) init <- c(init, orthAfter = 100)

    obj_fn <- function(pars) {
      para_sd <- rep(pars["paraBefore"], n_pairs)
      orth_sd <- rep(pars["orthBefore"], n_pairs)
      if (!is.na(para_break)) {
        para_sd[seq_len(n_pairs) > floor(para_break / 2)] <- pars["paraAfter"]
      }
      if (!is.na(orth_break)) {
        orth_sd[seq_len(n_pairs) > floor(orth_break / 2)] <- pars["orthAfter"]
      }
      -eval_ll(para_sd, orth_sd)
    }

    opt <- stats::optim(init, obj_fn, method = "L-BFGS-B",
                         lower = 0, upper = 1e10,
                         control = list(fnscale = 1))
    n_params <- 2 + (!is.na(para_break)) + (!is.na(orth_break))
    bic <- -2 * (-opt$value) + n_params * log(n)
    list(opt = opt, bic = bic)
  }

  res_none <- optimize_sds(NA, NA)
  best <- res_none
  best_pb <- NA
  best_ob <- NA
  for (pb in margin_breaks) {
    res <- optimize_sds(pb, NA)
    if (res$bic < best$bic) { best <- res; best_pb <- pb; best_ob <- NA }
  }
  for (ob in margin_breaks) {
    res <- optimize_sds(NA, ob)
    if (res$bic < best$bic) { best <- res; best_pb <- NA; best_ob <- ob }
  }
  if (!is.na(best_pb)) {
    for (ob in margin_breaks) {
      res <- optimize_sds(best_pb, ob)
      if (res$bic < best$bic) { best <- res; best_ob <- ob }
    }
  }
  if (!is.na(best_ob)) {
    for (pb in margin_breaks) {
      res <- optimize_sds(pb, best_ob)
      if (res$bic < best$bic) { best <- res; best_pb <- pb }
    }
  }
  opt <- best$opt
  result <- cbind(
    paraSd = rep(opt$par["paraBefore"], n),
    orthSd = rep(opt$par["orthBefore"], n)
  )
  if (!is.na(best_pb))
    result[seq_len(nrow(result)) >= best_pb, "paraSd"] <- opt$par["paraAfter"]
  if (!is.na(best_ob))
    result[seq_len(nrow(result)) >= best_ob, "orthSd"] <- opt$par["orthAfter"]
  min_brk <- min(c(margin_breaks), na.rm = TRUE)
  max_brk <- max(c(margin_breaks), na.rm = TRUE)
  result[seq_len(nrow(result)) < min_brk, ] <- NA
  result[seq_len(nrow(result)) >= max_brk, ] <- NA
  as.data.frame(result)
}

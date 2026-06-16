## Per-fix observation-error helpers shared across the package.
##
## location_error is a property of the location (a per-fix horizontal
## 1-sigma positional error, in metres).  Two consumers in the package:
##
## - mt_flag_outliers_bridge() injects the prior at the *anchors*
##   (i-1, i+1) to sharpen the bridge denominator while preserving
##   leverage immunity (the target fix's own sigma never enters its
##   own residual scale).  Anchor sigmas of NA fall back to the empirical
##   residual scale S_hat (no extra sharpening).
##
## - mt_dbbmm_variance() / mt_dbgb_variance() use the same per-fix
##   sigma as the prior anchor variance in the Brownian-bridge motion
##   model.  Here NA cannot fall back -- the per-fix variance enters
##   the diffusion integral directly -- so the variance entry points
##   take an explicit `location_error_na` strategy
##   ("median" / "mean" / "zero" / "approx").
##
## See mt_flag_outliers_bridge() and mt_dbbmm_variance() for the math;
## this file holds the input-resolver, the NA-imputation helper, and
## the Argos LC lookup table.


## Resolve user-supplied observation-error specification to a length-n
## vector of horizontal 1-sigma values in metres.
##
## location_error accepts:
##   NULL          - disabled, returns NULL
##   numeric(1)    - uniform sigma for all fixes
##   numeric(n)    - per-fix sigma already in metres
##   character(1)  - column name in x containing per-fix sigma in m,
##                   except for "auto" which probes known columns
##
## Returns NULL (disabled) or a length-n numeric vector with NA where
## sigma is unknown for a particular fix.  Downstream code treats NA
## as 0 (no obs-error contribution at that anchor).
##
## @keywords internal
.resolve_location_error <- function(location_error, x, n) {
  if (is.null(location_error)) return(NULL)

  if (is.numeric(location_error)) {
    if (length(location_error) == 1L) {
      if (is.na(location_error) || location_error < 0) {
        rlang::abort(
          "`location_error` scalar must be a non-negative finite number.",
          class = "move2utils_location_error_bad_scalar")
      }
      return(rep(as.numeric(location_error), n))
    }
    if (length(location_error) != n) {
      rlang::abort(
        sprintf("`location_error` numeric vector must have length nrow(x) = %d, got %d.",
                n, length(location_error)),
        class = "move2utils_location_error_bad_vector_length")
    }
    sig <- as.numeric(location_error)
    sig[!is.finite(sig) | sig < 0] <- NA_real_
    return(sig)
  }

  if (is.character(location_error) && length(location_error) == 1L) {
    if (location_error == "auto") {
      return(.resolve_location_error_auto(x, n))
    }
    col <- .find_col(location_error, names(x))
    if (is.na(col)) {
      rlang::abort(
        sprintf("`location_error` names a column (\"%s\") that is not in `x`.",
                location_error),
        class = "move2utils_location_error_missing_column")
    }
    sig <- suppressWarnings(as.numeric(x[[col]]))
    if (length(sig) != n) {
      rlang::abort(
        sprintf("Column \"%s\" has length %d, expected %d.",
                col, length(sig), n),
        class = "move2utils_location_error_bad_column_length")
    }
    sig[!is.finite(sig) | sig < 0] <- NA_real_
    return(sig)
  }

  rlang::abort(paste0(
    "`location_error` must be NULL, a numeric scalar, a numeric vector ",
    "of length nrow(x), or a character column name (or \"auto\")."),
    class = "move2utils_location_error_bad_type")
}


## Auto-detect a per-fix observation-error column.  Probes common
## Movebank columns in priority order and returns a length-n sigma
## vector in metres, or NULL with a message if nothing usable is
## found.
##
## Priority:
##   1. eobs_horizontal_accuracy_estimate (m, e-obs GPS)
##   2. argos_lc + .argos_lc_sigma() lookup
##
## gps_hdop / gps_dop need a user-known multiplier (the GPS receiver's
## URE) to convert to metres and so are not auto-resolved.  Users with
## those columns should construct sigma themselves and pass via column
## name.
##
## @keywords internal
.resolve_location_error_auto <- function(x, n) {
  cn <- names(x)

  hacc <- .find_col("eobs_horizontal_accuracy_estimate", cn)
  if (!is.na(hacc)) {
    sig <- suppressWarnings(as.numeric(x[[hacc]]))
    sig[!is.finite(sig) | sig < 0] <- NA_real_
    n_ok <- sum(!is.na(sig))
    message(sprintf(
      "  location_error = \"auto\": using `%s` (m); %d/%d fixes have a value.",
      hacc, n_ok, n))
    return(sig)
  }

  argos <- .find_col("argos_lc", cn)
  if (!is.na(argos)) {
    lc  <- as.character(x[[argos]])
    sig <- .argos_lc_sigma(lc)
    n_ok <- sum(!is.na(sig))
    message(sprintf(
      "  location_error = \"auto\": mapped `%s` to sigma via internal lookup; %d/%d fixes have a value.",
      argos, n_ok, n))
    return(sig)
  }

  message("  location_error = \"auto\": no recognised quality column found ",
          "(`eobs_horizontal_accuracy_estimate`, `argos_lc`); injection disabled.")
  NULL
}


## Impute NAs in a per-fix `location_error` vector using one of four
## strategies.  Caller (typically a variance _single fit) supplies the
## per-track slice; imputation is local to that slice.
##
## Strategies:
##   "median" - fill NAs with the per-track median of non-NA values.
##              Default for variance / UD entry points; conservative
##              central tendency, robust to a few extreme errors.
##   "mean"   - fill NAs with the per-track mean of non-NA values.
##   "zero"   - fill NAs with 0.  Fakes certainty -- only legitimate
##              when the analyst has explicitly decided to disregard
##              measurement error at NA fixes.
##   "approx" - linearly interpolate from neighbouring non-NA values
##              along the (positional) index, with edge extrapolation
##              via the nearest non-NA (rule = 2).
##
## Returns the vector with NAs filled.  If all values are NA, errors
## unless strategy = "zero" (where filling with 0 is well-defined).
##
## @keywords internal
.impute_loc_err_na <- function(loc_err, na_replace) {
  if (is.null(loc_err) || !any(is.na(loc_err))) return(loc_err)

  na_replace <- match.arg(na_replace,
                          c("median", "mean", "zero", "approx"))
  ok <- !is.na(loc_err)

  if (!any(ok) && na_replace != "zero") {
    rlang::abort(
      sprintf("location_error is entirely NA (n = %d). Cannot impute with na_replace = \"%s\". Pass na_replace = \"zero\" only if a uniform zero (no per-fix error contribution) is intentional.",
              length(loc_err), na_replace),
      class = "move2utils_location_error_all_na")
  }

  if (na_replace == "approx") {
    idx <- seq_along(loc_err)
    loc_err <- stats::approx(x = idx[ok], y = loc_err[ok],
                              xout = idx, rule = 2L)$y
    return(loc_err)
  }

  fill <- switch(na_replace,
                 median = stats::median(loc_err[ok]),
                 mean   = mean(loc_err[ok]),
                 zero   = 0)
  loc_err[!ok] <- fill
  loc_err
}


## Resolve a `location_error` specification *and* impute NAs in one
## call -- the entry point used by mt_dbbmm_variance and
## mt_dbgb_variance after they have sliced their per-track subset.
##
## Returns a length-`n` numeric vector with no NAs and no negatives.
## Scalar / vector inputs pass through .resolve_location_error first
## so column-name and "auto" forms are honoured uniformly.
##
## @keywords internal
.resolve_loc_err_for_variance <- function(location_error, x, n,
                                            na_replace = "median") {
  if (is.null(location_error)) return(NULL)

  if (is.numeric(location_error) && length(location_error) == 1L) {
    if (is.na(location_error) || location_error < 0) {
      rlang::abort(
        "`location_error` scalar must be a non-negative finite number.",
        class = "move2utils_location_error_bad_scalar")
    }
    return(rep(as.numeric(location_error), n))
  }

  sig <- .resolve_location_error(location_error, x, n)
  if (is.null(sig)) return(NULL)
  if (length(sig) == 1L) sig <- rep(sig, n)
  if (length(sig) != n) {
    rlang::abort(
      sprintf("resolved location_error has length %d, expected %d.",
              length(sig), n),
      class = "move2utils_location_error_resolved_length_mismatch")
  }

  if (any(is.na(sig))) {
    sig <- .impute_loc_err_na(sig, na_replace)
  }
  sig
}


## Map an Argos location class (a character vector with values in
## "3", "2", "1", "0", "A", "B", "Z") to a 1-sigma horizontal error in
## metres.
##
## Values follow the conventional CLS ladder cited in Vincent et al.
## (2002) and widely re-used in marine-mammal and seabird tracking
## studies (e.g. Boyd & Brightsmith 2013).  These are conservative
## central estimates; analysts with empirical calibrations specific to
## their tag, deployment, or region should pass their own column
## instead.
##
## Class Z is the rejected-fix flag and yields NA (no usable position).
##
## @keywords internal
.argos_lc_sigma <- function(lc) {
  tab <- c("3" = 250, "2" = 500, "1" = 1500, "0" = 5000,
           "A" = 5000, "B" = 10000, "Z" = NA_real_)
  lc <- toupper(as.character(lc))
  out <- rep(NA_real_, length(lc))
  ok <- !is.na(lc) & lc %in% names(tab)
  out[ok] <- tab[lc[ok]]
  out
}


## Return the column name in `cn` matching `underscore_name` either
## directly or with dashes (Movebank CSV convention vs. mt_read).
##
## @keywords internal
.find_col <- function(underscore_name, cn) {
  if (underscore_name %in% cn) return(underscore_name)
  dashed <- gsub("_", "-", underscore_name, fixed = TRUE)
  if (dashed %in% cn) return(dashed)
  NA_character_
}

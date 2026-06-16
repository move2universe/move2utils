## Central dispatcher for primitive-wrapper boilerplate
##
## Each public outlier-detection primitive (`mt_flag_outliers_bridge`,
## `mt_flag_outliers_detour`, ...) wraps a `.<primitive>_fn_core` raw-
## matrix entry point.  The wrapper has historically duplicated a fixed
## sequence of move2 -> .fn_core plumbing:
##
##   1. Multi-track dispatch (lapply over unique track ids)
##   2. AEQD projection if input is lon/lat
##   3. Edge case: n < n_min returns immediately with a warning
##   4. Extract `cc <- st_coordinates(x)`, `t_s <- mt_time(x)`,
##      `track_id <- mt_track_id(x)`
##   5. Input hygiene: non-finite coords / times; duplicate or out-of-
##      order timestamps (latter only when `hygiene_strict_time = TRUE`)
##   6. Build the initial active mask
##   7. Call `fn_core(cc, t_s, active_idx, ...)`
##   8. Lift active-indexed result vectors to full-length n; attach
##      result columns to x
##   9. CRS restore: re-attach result columns to the original (pre-
##      projection) move2 object so the caller's CRS is preserved
##
## `.clean_track_dispatch()` does all of this in one place.  Per-
## primitive wrappers become small shims that name their `.fn_core` +
## pass primitive-specific args.

#' Dispatch a `.fn_core` call from a `move2` object.
#'
#' Internal helper.  Handles the common plumbing each outlier-
#' detection primitive wrapper would otherwise duplicate.
#'
#' @param x A `move2` object.
#' @param fn_core The `.fn_core` function to call (e.g.
#'   `.bridge_fn_core`).
#' @param fn_core_args Named list of extra arguments to pass to
#'   `fn_core` *after* the standard positional ones (`cc`, optionally
#'   `t_s`, `active_idx`).
#' @param per_track_args Named list of arguments whose **numeric
#'   vector** values need to be sliced per track during multi-track
#'   dispatch (e.g. `location_error`).  Non-numeric or scalar values
#'   pass through unchanged.  After per-track slicing, these are
#'   merged into `fn_core_args` for the single-track call.
#' @param pool_args_by_track Named list, track id (as character) ->
#'   a named list of arguments to inject into `fn_core_args` for
#'   that track.  Used by `pool_by`: the wrapper fits one set of
#'   pooled parameters per group, then builds this map so every
#'   track in a group receives the same fitted parameters during
#'   dispatch.  NULL (default) means no pool injection.
#' @param need_time Logical.  If `TRUE` (default), extract a
#'   time-in-seconds vector and pass it to `fn_core` as the second
#'   positional argument.  If `FALSE`, omit time extraction (for
#'   geometric, time-insensitive primitives such as
#'   `mt_flag_outliers_detour`).
#' @param hygiene_strict_time Logical.  If `TRUE` and `need_time =
#'   TRUE`, stop with an informative error when consecutive duplicate
#'   or out-of-order timestamps are present.  If `FALSE`, skip the
#'   strict timestamp checks (still excludes rows with non-finite
#'   times from the active mask).
#' @param project_longlat Logical.  If `TRUE` (default) and the input
#'   is in geographic coordinates, project to a local AEQD before
#'   extracting `cc`; the result columns are then attached back to
#'   the original (geographic) `move2` object so the caller's CRS is
#'   preserved.  Set `FALSE` for primitives that handle geographic
#'   input directly (e.g. `mt_flag_outliers_detour` uses Haversine
#'   on the raw coordinates to avoid AEQD projection cost on
#'   multi-million-fix tracks).
#' @param n_min Integer.  Minimum row count below which the dispatcher
#'   emits a diagnostic (see `n_min_severity`) and calls `fn_core`
#'   anyway (which returns all-NA output via its built-in min-active
#'   guard).
#' @param n_min_severity One of `"warning"` (default; fires a `warning()`
#'   that the user sees regardless of `silent`) or `"say"` (a
#'   `message()` gated by `silent`).  Match the originating primitive's
#'   pre-refactor convention: e.g. `mt_flag_outliers_bridge` warned,
#'   `mt_flag_outliers_detour` used a quiet say-message.
#' @param primitive_label Character.  Short name for the primitive,
#'   used only in diagnostic messages (`"bridge"`, `"detour"`, ...).
#' @param silent Logical.  Suppress per-primitive narration.
#'
#' @return The augmented `move2` object with `fn_core`'s output
#'   columns attached (and any scalar elements of `fn_core`'s return
#'   list lifted to `attr(x, name)`).  CRS restored to the caller's
#'   original.
#'
#' @keywords internal
.clean_track_dispatch <- function(x,
                                  fn_core,
                                  fn_core_args        = list(),
                                  per_track_args      = list(),
                                  pool_args_by_track  = NULL,
                                  need_time           = TRUE,
                                  hygiene_strict_time = TRUE,
                                  project_longlat     = TRUE,
                                  n_min               = 10L,
                                  n_min_severity      = c("warning", "say"),
                                  primitive_label     = "primitive",
                                  silent              = FALSE) {
  n_min_severity <- match.arg(n_min_severity)

  say <- function(...) if (!silent) message(...)

  ## Empty geometry is rejected up front (no outlier in a fix with no
  ## location); checked here before any coordinate extraction or track split.
  .reject_empty_geometry(x, sprintf("%s scoring", primitive_label))

  ## ---- 1. Multi-track dispatch ----
  ids <- move2::mt_track_id(x)
  unique_ids <- unique(ids)
  if (length(unique_ids) > 1L) {
    say(sprintf("Processing %d individuals separately...",
                length(unique_ids)))
    results <- lapply(unique_ids, function(id) {
      idx <- which(ids == id)
      say(sprintf("--- %s (%d locations) ---", as.character(id),
                  length(idx)))
      pta_id <- lapply(per_track_args, function(v) {
        if (is.atomic(v) && !is.null(v) && length(v) > 1L) v[idx] else v
      })
      ## Inject pooled args for this track if pool_args_by_track is set.
      ## Each track in a pool group receives the same fitted parameters.
      pool_id <- if (!is.null(pool_args_by_track)) {
        pool_args_by_track[[as.character(id)]]
      } else list()
      if (is.null(pool_id)) pool_id <- list()
      .clean_track_dispatch(x[idx, , drop = FALSE],
                              fn_core             = fn_core,
                              fn_core_args        = c(fn_core_args, pta_id, pool_id),
                              per_track_args      = list(),
                              pool_args_by_track  = NULL,
                              need_time           = need_time,
                              hygiene_strict_time = hygiene_strict_time,
                              project_longlat     = project_longlat,
                              n_min               = n_min,
                              n_min_severity      = n_min_severity,
                              primitive_label     = primitive_label,
                              silent              = silent)
    })
    return(do.call(rbind, results))
  }

  ## ---- 1b. Single-track call with pool_args_by_track ----
  ## The wrapper may have called us on a single-track input but with
  ## a pool_args_by_track map (from a parent multi-track recursion
  ## that has been collapsed by an upstream caller).  In that case,
  ## look up this track's id and inject.
  if (!is.null(pool_args_by_track)) {
    sole_id <- as.character(unique_ids[[1]])
    pool_id <- pool_args_by_track[[sole_id]]
    if (!is.null(pool_id)) {
      fn_core_args <- c(fn_core_args, pool_id)
    }
  }

  ## ---- 2. Canonicalise to a local AEQD (when project_longlat) ----
  ## Project to the per-track canonical AEQD regardless of the input CRS, so
  ## the geometry is computed in metres and the result is invariant to the
  ## CRS the caller supplied.  (Previously only longlat input was projected;
  ## a caller-supplied projection was used as-is, which leaked their CRS into
  ## the result.)  Detectors whose math is already scale/CRS-invariant
  ## (detour ratio, speed cap via metric distances) pass project_longlat =
  ## FALSE and stay in the input CRS to avoid the projection cost.
  orig_was_longlat <- isTRUE(sf::st_is_longlat(x))
  orig_x <- x
  did_project <- isTRUE(project_longlat)
  if (did_project) {
    if (orig_was_longlat) {
      say(sprintf(paste0(
        "Input is in longitude/latitude.  Auto-projecting to a local ",
        "AEQD for Euclidean %s math; output is returned in the ",
        "original CRS."), primitive_label))
    }
    x <- .to_canonical_aeqd(x)
    ## geometry is now metric; tell the core not to use longlat formulae
    if ("was_longlat" %in% names(fn_core_args)) fn_core_args$was_longlat <- FALSE
  }

  n <- nrow(x)

  ## ---- 3. Edge case: n < n_min ----
  ## .fn_core handles the actual empty-output case via its internal
  ## min-active guard; we just emit the diagnostic.  Severity per
  ## the originating primitive's convention (some warn, some
  ## silently-message).
  if (n < n_min) {
    msg <- sprintf("Too few locations for %s (<%d). Returning without flags.",
                   primitive_label, n_min)
    if (n_min_severity == "warning") {
      rlang::warn(msg,
                  class = "move2utils_primitive_too_few_locations")
    } else {
      say(msg)
    }
  }

  ## ---- 4. Extract cc, t_s, track_id (already a single id here) ----
  cc_all <- sf::st_coordinates(x)
  t_all  <- if (need_time) as.numeric(move2::mt_time(x), units = "secs") else NULL

  ## ---- 5. Input hygiene ----
  bad_coord <- !is.finite(cc_all[, 1]) | !is.finite(cc_all[, 2])
  bad_time  <- if (need_time) !is.finite(t_all) else rep(FALSE, n)
  bad_row   <- bad_coord | bad_time

  if (need_time && hygiene_strict_time) {
    dt_all   <- diff(t_all)
    n_dup_t  <- sum(dt_all == 0, na.rm = TRUE)
    n_neg_dt <- sum(dt_all < 0,  na.rm = TRUE)
    if (n_dup_t > 0 || n_neg_dt > 0) {
      msgs <- character(0)
      if (n_dup_t > 0) {
        msgs <- c(msgs, sprintf(
          "  - %d consecutive pair(s) with duplicate timestamps. Run `move2::mt_filter_unique(x, criterion = \"first\")` (or \"sample\", \"fix_type\", etc.) before calling this function.",
          n_dup_t))
      }
      if (n_neg_dt > 0) {
        msgs <- c(msgs, sprintf(
          "  - %d consecutive pair(s) with negative Delta-t (track is not time-sorted). Sort by `mt_time(x)` within each track before calling.",
          n_neg_dt))
      }
      ## Signal under the shared precondition classes (so handlers keyed on
      ## move2utils_input_duplicate_timestamps / _unsorted_timestamps fire
      ## here too), keeping the umbrella class as a secondary identifier.
      anomaly_classes <- c(
        if (n_dup_t > 0)  "move2utils_input_duplicate_timestamps",
        if (n_neg_dt > 0) "move2utils_input_unsorted_timestamps",
        "move2utils_temporal_anomaly")
      rlang::abort(
        paste0(
          sprintf("Input has temporal anomalies that prevent %s scoring:\n",
                  primitive_label),
          paste(msgs, collapse = "\n"),
          "\nResolve these upstream, then re-run."),
        class = anomaly_classes)
    }
  }

  ## Non-finite coordinates (NA geometry) are rejected, like empty geometry:
  ## there is no outlier to identify in a fix with no usable location.
  if (any(bad_coord)) {
    rlang::abort(
      sprintf(paste0(
        "%d fix(es) have non-finite coordinates. Remove them before %s ",
        "scoring, e.g. `mt_filter_gps_quality(x)`. There is no outlier to ",
        "identify in a fix with no usable location."),
        sum(bad_coord), primitive_label),
      class = "move2utils_input_nonfinite_coords")
  }
  ## Non-finite timestamps cannot be scored (no time lag); these are
  ## excluded from scoring with a note rather than rejected.
  if (need_time && any(bad_time)) {
    say(sprintf(paste0(
      "Input-hygiene note: %d fix(es) with non-finite timestamps are ",
      "excluded from %s scoring. To preprocess upstream, use ",
      "`mt_filter_gps_quality(x)`."), sum(bad_time), primitive_label))
  }

  ## ---- 6. Initial active mask ----
  active_idx_init <- which(!bad_row)

  ## ---- 7. Call .fn_core ----
  ## Single-track path: per_track_args carry full-length values that
  ## need no slicing (the input was single-track to begin with), so
  ## fold them straight into the fn_core call.
  base_args <- if (need_time) {
    list(cc_all, t_all, active_idx_init)
  } else {
    list(cc_all, active_idx_init)
  }
  res <- do.call(fn_core, c(base_args, fn_core_args, per_track_args))

  ## ---- 8. Lift active-indexed result to full-length n -----
  ## Convention: scalar elements (length 1) are lifted to attributes;
  ## non-scalar elements are lifted as full-length-n columns of the
  ## same type, with NA / FALSE outside `active_idx_init`.
  for (nm in names(res)) {
    val <- res[[nm]]
    if (is.null(val)) next
    if (length(val) == 1L && !is.atomic(val[[1]])) next  # unlikely
    if (length(val) == 1L) {
      attr(x, nm) <- val
      next
    }
    if (length(val) != length(active_idx_init)) {
      rlang::abort(
        sprintf(paste0(
          "Internal error: %s_fn_core's `%s` has length %d but ",
          "active_idx has length %d.  This is a bug in the primitive's ",
          ".fn_core; the dispatcher expects active-indexed return ",
          "vectors."), primitive_label, nm, length(val),
          length(active_idx_init)),
        class = "move2utils_internal_fn_core_bad_active_length")
    }
    fill <- if (is.logical(val)) {
              if (nm == "is_outlier") FALSE else NA
            } else if (is.integer(val)) {
              NA_integer_
            } else if (is.character(val)) {
              NA_character_
            } else {
              NA_real_
            }
    full <- rep(fill, n)
    full[active_idx_init] <- val
    x[[nm]] <- full
  }

  ## ---- 9. CRS restore (only when we actually projected) ----
  if (did_project) {
    new_cols  <- setdiff(names(x), names(orig_x))
    for (col in new_cols) orig_x[[col]] <- x[[col]]
    new_attrs <- setdiff(names(attributes(x)), names(attributes(orig_x)))
    for (a in new_attrs) attr(orig_x, a) <- attr(x, a)
    x <- orig_x
  }

  ## ---- 10. Per-track summary ----
  if (!is.null(x$is_outlier)) {
    say(sprintf("=== %d outliers (%.2f%% of %d) ===",
                sum(x$is_outlier), 100 * sum(x$is_outlier) / n, n))
  }

  x
}

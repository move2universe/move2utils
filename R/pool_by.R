## Shared infrastructure for the `pool_by` argument used by every
## threshold-fitting outlier primitive (mt_flag_outliers,
## mt_flag_outliers_bridge, mt_flag_outliers_detour, mt_flag_speed_cap)
## and the mt_clean_track orchestrator.
##
## `pool_by = NULL` (default everywhere) preserves the legacy
## per-track-id behaviour: each track is its own statistical universe
## and is processed independently in `.clean_track_dispatch`'s
## per-track lapply.
##
## `pool_by = "<column>"` names a column in `mt_track_data(x)` and
## groups tracks sharing that column's value.  The named column is
## used both as the source of the threshold-fitting distribution
## (events are unioned across tracks in the same group) and as the
## scope of the post-cascade flag union.  Under this single-column
## form, fit and union coincide.
##
## `pool_by = c("<outer>", "<inner>")` separates the two roles:
##   * `outer`  names the column whose value defines the *fit set*:
##                threshold-fitting primitives draw their reference
##                distribution from the union of events whose tracks
##                share this column's value.
##   * `inner`  names the column whose value defines the *operating
##                unit*: the post-cascade flag union acts within tracks
##                sharing this column's value (and the orchestrator's
##                post-cascade sweep iterates over these groups).
## The two-element form requires `inner` to nest strictly within
## `outer`: every distinct `inner` value must map to exactly one
## `outer` value.  Length-1 input behaves identically to
## `c(outer = X, inner = X)`.
##
## Cap at length 2: the API has exactly two semantic roles
## (distribution source vs. operating unit).  A deeper hierarchy
## (e.g. species/population/individual/tag) would only earn its keep
## under hierarchical / partial-pooling threshold estimation, which
## the cascade does not perform.  When users have many nested levels,
## they pick the *pair* of columns that captures their trust claim
## (which level's distribution to fit from) and their operating unit
## (which level to union over).
##
## Per-primitive pool-fit logic lives in each primitive's wrapper,
## since "what to fit" differs (bridge_residual quantiles, joint
## step/turn/autodiff densities, speed quantiles, etc.).  The
## wrappers consume the helper's outer / inner maps and either build
## a `pool_args_by_track` named list (integrated path, used by
## `mt_flag_outliers`) or pass a `pool_step_fn(fit_idx, apply_idx)`
## closure to `.apply_pool_union` (post-hoc path, used by the other
## three primitives).


#' Resolve `pool_by` into outer and inner track-id -> group-id maps.
#'
#' Internal helper.  Validates the named column(s) exist in the
#' move2 object's track-level metadata, returns a two-element list
#' of named character maps from track id (as character) to group id
#' (as character).
#'
#' Tracks whose `pool_by` column value is NA fall back to their own
#' one-track group (group id = track id) with a one-line warning.
#' This is applied per map (outer and inner each handle NA
#' independently).
#'
#' @param x A move2 object.
#' @param pool_by Character vector of length 1 or 2 naming columns
#'   in `mt_track_data(x)`.  Length 1: same column is used as both
#'   outer (fit set) and inner (operating unit).  Length 2:
#'   `c(outer, inner)` -- `outer` defines the threshold-fit
#'   distribution source, `inner` defines the flag-union /
#'   post-cascade sweep unit.  `inner` must nest strictly in `outer`:
#'   every distinct `inner` value must map to exactly one `outer`
#'   value.  Anything longer is rejected with the
#'   hierarchical-estimation rationale.  Must be non-NULL; callers
#'   gate the NULL case before invoking this helper.
#' @param silent Logical.  Suppress the NA-fallback warning.
#'
#' @return Named list with two elements:
#'   \describe{
#'     \item{`outer`}{Named character vector mapping track id ->
#'       outer-group id, used as the fit-set selector.}
#'     \item{`inner`}{Named character vector mapping track id ->
#'       inner-group id, used as the union / operating-unit selector.}
#'   }
#'   Both maps share the same set of names (all unique track ids of
#'   `x`, as character).  Under length-1 input, `outer` and `inner`
#'   are identical.
#'
#' @keywords internal
.resolve_pool_groups <- function(x, pool_by, silent = FALSE) {
  if (is.null(pool_by)) {
    rlang::abort(paste0(
      "Internal error: .resolve_pool_groups called with pool_by = NULL. ",
      "Callers must gate the NULL case before invoking."),
      class = "move2utils_internal_pool_by_null_in_resolve")
  }
  if (!is.character(pool_by) || !length(pool_by) %in% c(1L, 2L)) {
    rlang::abort(
      sprintf(paste0(
        "`pool_by` must be NULL, or a character vector of length 1 ",
        "(single column: same column used as both fit set and operating ",
        "unit) or length 2 (`c(outer, inner)`: outer = distribution ",
        "source, inner = operating unit). Deeper hierarchies (>2 levels) ",
        "would require hierarchical threshold estimation, which the ",
        "cascade does not perform. If your data has more nested levels ",
        "(e.g. species/population/individual/tag), pick the two columns ",
        "that match your scientific claim: which level's distribution do ",
        "you trust, and which level should the flag union respect? ",
        "Got: %s"),
        deparse(pool_by)),
      class = "move2utils_pool_by_bad_type")
  }
  if (length(pool_by) == 2L && pool_by[[1L]] == pool_by[[2L]]) {
    rlang::abort(
      sprintf(paste0(
        "`pool_by` of length 2 must name two *distinct* columns ",
        "(`c(outer, inner)`). Got identical values: %s. ",
        "Use length-1 form instead."),
        deparse(pool_by[[1L]])),
      class = "move2utils_pool_by_duplicate_columns")
  }

  td <- move2::mt_track_data(x)
  missing <- setdiff(pool_by, names(td))
  if (length(missing)) {
    rlang::abort(
      paste0(
        sprintf("pool_by names column(s) not in `mt_track_data(x)`: %s. ",
                paste(shQuote(missing), collapse = ", ")),
        sprintf("Available columns: %s",
                paste(shQuote(names(td)), collapse = ", "))),
      class = "move2utils_pool_by_missing_columns")
  }

  ## mt_track_data row order matches mt_n_tracks() — one row per
  ## unique track_id, in some package-defined order.  Extract the
  ## track-id column explicitly so we can match by name rather than
  ## relying on row position.
  tid_col <- move2::mt_track_id_column(x)
  if (!(tid_col %in% names(td))) {
    ids_unique <- as.character(unique(move2::mt_track_id(x)))
    if (length(ids_unique) != nrow(td)) {
      rlang::abort(
        paste0(
          sprintf("Cannot resolve pool_by: the track-id column \"%s\" ",
                  tid_col),
          "is not in mt_track_data(x) and the unique-id count ",
          sprintf("(%d) does not match mt_track_data() row count (%d). ",
                  length(ids_unique), nrow(td)),
          "Please file a bug report."),
        class = "move2utils_internal_pool_by_id_count_mismatch")
    }
    ids_for_join <- ids_unique
  } else {
    ids_for_join <- as.character(td[[tid_col]])
  }

  resolve_one <- function(col, role) {
    raw <- as.character(td[[col]])
    na_mask <- is.na(raw) | !nzchar(raw)
    out <- raw
    if (any(na_mask)) {
      out[na_mask] <- ids_for_join[na_mask]   # NA -> own one-track group
      if (!silent) {
        rlang::warn(
          paste0(
            sprintf("%d track(s) have NA in pool_by %s column \"%s\"; processing each independently. ",
                    sum(na_mask), role, col),
            "To pool these, ensure the column is populated upstream."),
          class = "move2utils_pool_by_na_column")
      }
    }
    names(out) <- ids_for_join
    out
  }

  outer_col <- pool_by[[1L]]
  inner_col <- if (length(pool_by) == 2L) pool_by[[2L]] else pool_by[[1L]]

  outer_map <- resolve_one(outer_col, role = if (length(pool_by) == 2L) "outer" else "")
  inner_map <- if (length(pool_by) == 2L) {
    resolve_one(inner_col, role = "inner")
  } else {
    outer_map
  }

  ## Strict-nesting validation (length-2 only): every distinct inner
  ## value must map to a single outer value.  Single-level form is
  ## trivially self-nested (outer == inner).
  if (length(pool_by) == 2L) {
    inner_to_outer <- split(unname(outer_map), inner_map)
    bad <- inner_to_outer[vapply(inner_to_outer,
                                 function(v) length(unique(v)) > 1L,
                                 logical(1))]
    if (length(bad)) {
      sample_bad <- names(bad)[1L]
      rlang::abort(
        paste0(
          sprintf("`pool_by = c(\"%s\", \"%s\")` requires inner to nest in outer, ",
                  outer_col, inner_col),
          sprintf("but inner value \"%s\" spans %d distinct outer values: %s. ",
                  sample_bad, length(unique(bad[[1L]])),
                  paste(shQuote(unique(bad[[1L]])), collapse = ", ")),
          "Every distinct inner value must map to exactly one outer value. ",
          "Either fix the metadata so the hierarchy is consistent, or ",
          "choose a different (outer, inner) pair."),
        class = "move2utils_pool_by_not_nested")
    }
  }

  list(outer = outer_map, inner = inner_map)
}


#' Apply a pool-fit threshold + union flags into an existing output object.
#'
#' Internal helper shared by every post-hoc threshold-fitting primitive
#' (`mt_flag_outliers_detour`, `mt_flag_speed_cap`,
#' `mt_flag_outliers_bridge`).  The primitive's wrapper runs the
#' normal per-track dispatch first (producing per-track flags with
#' track-local thresholds), then calls this helper to layer a
#' pool-fitted threshold on top -- additive only, never un-flags
#' anything the per-track pass caught.
#'
#' Iterates over INNER groups (the operating unit).  For each inner
#' group, the OUTER group it nests in supplies the threshold-fit
#' event set; the inner group itself supplies the events the new
#' threshold is evaluated against and unioned into.  Under length-1
#' `pool_by`, outer == inner and both index sets coincide --
#' identical to the legacy single-level behaviour.
#'
#' Single-track-in-group equivalence: if an inner group contains one
#' track AND its outer group also contains only that track, the
#' pool-fit threshold is computed from that track's own diagnostic,
#' which is the same data the per-track threshold was computed from.
#' The pool union is then by construction byte-identical to the
#' per-track flags (assuming the `pool_step_fn` reproduces the per-
#' track gates exactly).  Pool_by never regresses a single-track
#' input.
#'
#' @param out The output `move2` object from per-track dispatch.
#'   Must already carry the diagnostic columns the `pool_step_fn`
#'   closure will read.
#' @param x The original input `move2` (for `mt_track_id` and
#'   `mt_track_data`).
#' @param pool_by Character vector length 1 or 2 (already validated
#'   by caller).
#' @param pool_step_fn Function `function(fit_idx, apply_idx) ->
#'   logical(length(apply_idx))`.  Given integer event-row indices
#'   of the OUTER group (`fit_idx`, the fit set) and INNER group
#'   (`apply_idx`, the union target), fit the pool threshold from
#'   `fit_idx`'s data and return a logical vector marking which
#'   `apply_idx` events the pool threshold would flag.  Under
#'   length-1 `pool_by`, `fit_idx` and `apply_idx` are the same
#'   index set.  The closure is responsible for primitive-specific
#'   gating (`legs_ok`, mode-position guard, etc.) so the union
#'   respects all gates the per-track flagging would have applied.
#' @param flag_cols Character vector.  Names of logical columns
#'   in `out` to union the pool flags into (e.g.
#'   `c("is_outlier", "flagged_by_detour")`).  Missing columns
#'   are silently skipped.
#' @param silent Logical.  Suppress the NA-fallback warning from
#'   `.resolve_pool_groups`.
#'
#' @return `out` with `flag_cols` unioned in place.
#'
#' @keywords internal
.apply_pool_union <- function(out, x, pool_by, pool_step_fn,
                              flag_cols, silent = FALSE) {
  maps <- .resolve_pool_groups(x, pool_by, silent = silent)
  outer_map <- maps$outer
  inner_map <- maps$inner

  ## Use mt_track_id on `x`, not `out` -- the dispatcher may rbind
  ## tracks in a different row order than `x`, so we have to match
  ## by track id rather than by row position.  The lift preserves
  ## the original row order within each track, but the inter-track
  ## order is dispatcher-controlled.
  ids_event <- as.character(move2::mt_track_id(out))

  ## pool_added records, per event row, whether the pool path added
  ## a flag that was NOT already in `is_outlier` from the per-track
  ## dispatch.  Used by the orchestrator's post-cascade sweep to
  ## identify pool-contributed flags without re-running the primitive
  ## without pool_by.  Length = nrow(out).
  pool_added <- logical(nrow(out))
  prior_outlier <- if ("is_outlier" %in% names(out)) {
    v <- out$is_outlier
    v[is.na(v)] <- FALSE
    v
  } else {
    logical(nrow(out))
  }

  for (g_inner in unique(inner_map)) {
    inner_tracks <- names(inner_map)[inner_map == g_inner]
    if (!length(inner_tracks)) next
    apply_idx <- which(ids_event %in% inner_tracks)
    if (!length(apply_idx)) next

    ## Strict nesting guarantees every inner-track shares one outer
    ## value, so taking the outer of any inner track is sufficient.
    g_outer <- unname(outer_map[inner_tracks[[1L]]])
    outer_tracks <- names(outer_map)[outer_map == g_outer]
    fit_idx <- which(ids_event %in% outer_tracks)
    if (!length(fit_idx)) next

    new_flag <- pool_step_fn(fit_idx, apply_idx)
    if (length(new_flag) != length(apply_idx)) {
      rlang::abort(
        sprintf("Internal error in pool_step_fn: returned %d flags for %d events.",
                length(new_flag), length(apply_idx)),
        class = "move2utils_internal_pool_step_fn_bad_length")
    }
    new_flag[is.na(new_flag)] <- FALSE
    pool_added[apply_idx] <- new_flag & !prior_outlier[apply_idx]
    for (col in flag_cols) {
      if (col %in% names(out)) {
        cur <- out[[col]][apply_idx]
        cur[is.na(cur)] <- FALSE
        out[[col]][apply_idx] <- cur | new_flag
      }
    }
  }
  attr(out, "pool_added") <- pool_added
  out
}


#' Validate state-vocabulary consistency within each inner-pool group.
#'
#' Internal helper used by `mt_clean_track` when both `pool_by` and
#' `state =` are supplied.  For each *inner* pool group (the
#' operating unit) containing more than one track, this checks that
#' every pairwise non-NA state vocabulary intersection within the
#' group is non-empty.  Outer-group consistency is not enforced --
#' the outer group is a fit source, not a union target.
#'
#' "Missing states" (a state in track A that never appears in
#' track B) is fine -- only entirely disjoint vocabularies count
#' as a mismatch.
#'
#' @param state_by_track Named list, track id (as character) -> the
#'   vector of state labels observed in that track (after NA strip).
#' @param inner_map Named character vector from `.resolve_pool_groups()$inner`.
#'
#' @return Invisibly TRUE on success; raises an informative error
#'   on disjoint vocabulary.
#'
#' @keywords internal
.validate_state_pool_vocab <- function(state_by_track, inner_map) {
  unique_groups <- unique(inner_map)
  for (g in unique_groups) {
    g_tracks <- names(inner_map)[inner_map == g]
    if (length(g_tracks) <= 1L) next   # nothing to compare
    vocabs <- lapply(g_tracks, function(t) {
      s <- state_by_track[[t]]
      s <- s[!is.na(s)]
      if (length(s)) unique(s) else character(0)
    })
    names(vocabs) <- g_tracks
    ## Drop empty vocabularies (tracks where state is all-NA) from the
    ## comparison -- they can't mismatch with anything.
    keep <- lengths(vocabs) > 0L
    vocabs <- vocabs[keep]
    g_tracks_eff <- names(vocabs)
    if (length(g_tracks_eff) <= 1L) next

    for (i in seq_along(g_tracks_eff)) {
      for (j in seq_len(i - 1L)) {
        inter <- intersect(vocabs[[i]], vocabs[[j]])
        if (!length(inter)) {
          rlang::abort(
            sprintf(paste0(
              "State vocabulary mismatch in pool unit \"%s\" between ",
              "tracks \"%s\" (states: %s) and \"%s\" (states: %s). ",
              "These vocabularies are entirely disjoint, which suggests ",
              "different state labelling conventions across tracks in ",
              "the same pool unit.  Resolve upstream (use the same state ",
              "labels across all tracks of one pool unit), or split the ",
              "tracks into different pool units, or run without ",
              "`pool_by`."),
              g,
              g_tracks_eff[[i]], paste(vocabs[[i]], collapse = ", "),
              g_tracks_eff[[j]], paste(vocabs[[j]], collapse = ", ")),
            class = "move2utils_pool_by_state_vocab_mismatch")
        }
      }
    }
  }
  invisible(TRUE)
}

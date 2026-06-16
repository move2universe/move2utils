#' Thin a track by along-track distance
#'
#' Select a subset of fixes spaced by cumulative along-track distance,
#' returning the `move2` object whole with a logical column
#' `thin_selected` marking retained fixes. Two methods are offered:
#'
#' * `method = "step"` (default) walks the track from the first fix and
#'   retains a fix each time the cumulative distance **since the last
#'   retained fix** reaches `distance`. Every retained pair is therefore
#'   at least `distance` of travel apart -- the intuitive "give me one
#'   fix per `distance` of travel" behaviour.
#' * `method = "interval"` is the distance-space twin of
#'   [mt_thin_time()]: it finds the *largest* subset such that every
#'   successive retained pair is within `[distance - tolerance,
#'   distance + tolerance]` of along-track travel. This is the faithful
#'   `move2` analogue of `move::thinDistanceAlongTrack()`, which is
#'   itself the distance counterpart of `move::thinTrackTime()`.
#'
#' @details
#' Per-segment distance comes from `move2::mt_distance(units = "m")`, so
#' `distance` (and `tolerance`) are always in metres regardless of the
#' input CRS's linear unit. Geographic (lon/lat) data are thinned
#' against great-circle distance and projected data against planar
#' distance, so the two can differ slightly at the same metre threshold;
#' project to an equal-area metric CRS first if you need the result to be
#' identical across coordinate systems. Multi-track input is thinned
#' independently per individual.
#'
#' For `method = "step"`, the first fix of each track is always retained,
#' and -- because the threshold is crossed mid-segment and we do not
#' interpolate -- the realised distance between successive retained fixes
#' can exceed `distance` (it never falls below it). Use
#' [move2::mt_interpolate()] first if you need exact-distance spacing.
#'
#' For `method = "interval"`, the track is pre-split at any single step
#' longer than `distance + tolerance` (thinning never bridges such a
#' step), and the longest tolerance-admissible chain is selected by the
#' same dynamic-programming sweep used by [mt_thin_time()]. The first fix
#' is not guaranteed to be retained, and no retained pair is closer than
#' `distance - tolerance`.
#'
#' @param x A `move2` object.
#' @param distance The target along-track distance step, in metres.
#'   Accepts a numeric scalar or a `units` object convertible to metres
#'   (e.g. `units::set_units(0.5, "km")`).
#' @param tolerance Half-width of the acceptance window around
#'   `distance`, in metres. Used only by `method = "interval"`; ignored
#'   (with a warning if set) by `method = "step"`. Defaults to
#'   `distance / 10`.
#' @param criterion Tie-breaking rule for `method = "interval"` when
#'   multiple chains of equal length exist: `"closest"` (minimum total
#'   deviation from the target, default) or `"first"` (earliest-starting
#'   chain). Ignored by `method = "step"`.
#' @param method Thinning rule: `"step"` (default, greedy minimum-spacing)
#'   or `"interval"` (tolerance-constrained longest chain).
#' @param remove Logical. If `TRUE`, return only the retained fixes.
#'   Default `FALSE` -- the full object is returned with a logical
#'   `thin_selected` column.
#'
#' @return A `move2` object. With `remove = FALSE`, the input is
#'   returned unchanged except for a new logical column
#'   `thin_selected`. With `remove = TRUE`, only rows with
#'   `thin_selected == TRUE` are returned.
#'
#' @examples
#' \donttest{
#' library(move2)
#' fishers <- mt_read(mt_example())
#' fishers <- fishers[!sf::st_is_empty(fishers), ]
#' leroy <- fishers[mt_track_id(fishers) == "M4", ][seq_len(500), ]
#'
#' ## one fix per ~300 m of travel (minimum spacing guaranteed)
#' out <- mt_thin_distance(leroy, distance = 300)
#' table(out$thin_selected)
#'
#' ## tolerance-constrained: retained pairs ~300 m +/- 30 m apart
#' out2 <- mt_thin_distance(leroy, distance = 300, tolerance = 30,
#'                          method = "interval")
#' table(out2$thin_selected)
#' }
#'
#' @seealso [mt_thin_time()] for time-based thinning (same DP machinery);
#'   [move2::mt_interpolate()] to first regularise the track in time or
#'   along a target line; [move2::mt_distance()] for the underlying
#'   distance.
#'
#' @export
mt_thin_distance <- function(x,
                             distance,
                             tolerance = NULL,
                             criterion = c("closest", "first"),
                             method    = c("step", "interval"),
                             remove    = FALSE) {
  stopifnot(inherits(x, "move2"))
  method    <- match.arg(method)
  criterion <- match.arg(criterion)

  if (inherits(distance, "units")) {
    distance <- as.numeric(units::set_units(distance, "m"))
  }
  stopifnot(is.numeric(distance), length(distance) == 1, distance > 0)

  ## ---- resolve tolerance (interval mode only) --------------------
  tol <- NULL
  if (method == "interval") {
    if (is.null(tolerance)) {
      tol <- distance / 10
    } else {
      if (inherits(tolerance, "units")) {
        tolerance <- as.numeric(units::set_units(tolerance, "m"))
      }
      stopifnot(is.numeric(tolerance), length(tolerance) == 1)
      if (tolerance < 0) {
        rlang::abort("`tolerance` must be non-negative.",
                     class = "move2utils_mt_thin_distance_bad_tolerance")
      }
      if (tolerance > distance) {
        rlang::abort("`tolerance` must not exceed `distance`.",
                     class = "move2utils_mt_thin_distance_tolerance_exceeds_distance")
      }
      tol <- tolerance
    }
  } else if (!is.null(tolerance)) {
    rlang::warn(
      "`tolerance` is ignored when `method = \"step\"`; it applies only to `method = \"interval\"`.",
      class = "move2utils_mt_thin_distance_tolerance_ignored")
  }

  n <- nrow(x)
  if (n == 0) {
    x$thin_selected <- logical(0)
    return(if (remove) x else x)
  }

  track_ids <- as.character(move2::mt_track_id(x))
  selected  <- logical(n)

  for (id in unique(track_ids)) {
    idx <- which(track_ids == id)
    if (length(idx) == 0L) next
    ## units = "m" so distances match the metres-normalised `distance`
    ## threshold regardless of the track's CRS (geodesic for lon/lat,
    ## converted for non-metre projected CRSs); as.numeric() alone would
    ## strip the unit without converting -- a CRS-unit leak.
    d <- as.numeric(move2::mt_distance(x[idx, ], units = "m"))
    d[is.na(d)] <- 0
    ## cumulative distance TO each fix (monotone non-decreasing)
    cumd <- cumsum(c(0, utils::head(d, -1L)))
    selected[idx] <- if (method == "step") {
      .thin_dist_step(cumd, distance)
    } else {
      ## cumulative distance is monotone, so the timestamp DP of
      ## mt_thin_time() applies verbatim with distance playing the role
      ## of the target interval.
      .thin_dp_track(cumd, distance, tol, criterion)
    }
  }

  x$thin_selected <- selected
  if (remove) x[selected, ] else x
}

## ------------------------------------------------------------------
## Greedy minimum-spacing walk: keep the first fix, then keep a fix each
## time cumulative travel since the last retained fix reaches `distance`.
.thin_dist_step <- function(cumd, distance) {
  m <- length(cumd)
  if (m == 0L) return(logical(0))
  keep <- logical(m)
  keep[1L] <- TRUE
  last <- cumd[1L]
  for (i in seq_len(m)[-1L]) {
    if (cumd[i] - last >= distance) {
      keep[i] <- TRUE
      last <- cumd[i]
    }
  }
  keep
}

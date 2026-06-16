#' Thin a track to a target time interval (tolerance-constrained)
#'
#' Select a subset of fixes such that every successive retained pair has
#' a time lag within `[interval - tolerance, interval + tolerance]`, and
#' such that the number of retained fixes is maximised. The original
#' `move2` object is returned whole with a logical column
#' `thin_selected` marking retained fixes.
#'
#' This is the move2 analogue of `move::thinTrackTime()`. The legacy
#' implementation used memoised recursion on 100-fix chunks; this one
#' uses a linear-time dynamic-programming sweep, so it scales cleanly
#' to tracks with \eqn{10^5}{1e5}+ fixes.
#'
#' @details
#' The problem is to find the longest path in a directed acyclic graph
#' whose nodes are fixes and whose edges `i -> j` exist iff
#' `t[j] - t[i]` is within `[interval - tolerance, interval + tolerance]`.
#' On sorted timestamps it is solved with the recursion
#' \deqn{f(i) = 1 + \max_{j < i, \; t_i - t_j \in [d - \tau, d + \tau]} f(j)}
#' followed by a single backtrace through the predecessor array.
#'
#' Tracks are pre-split at any gap longer than `interval + tolerance`;
#' thinning never bridges such a gap.
#'
#' When several chains are equally long, `criterion = "closest"` (the
#' default) picks the chain whose summed absolute deviation from the
#' target interval is smallest; `criterion = "first"` picks the chain
#' that starts earliest. The `"all"` mode of the legacy function is not
#' reproduced — it can explode combinatorially and is rarely the
#' operationally useful choice. If you need exact parity with
#' `move::thinTrackTime(..., criteria = "all")`, call the legacy
#' function on a `move2::to_move()` conversion.
#'
#' @param x A `move2` object.
#' @param interval Target interval between successive retained fixes.
#'   A `difftime`, a `units` object with a time dimension (e.g.
#'   `units::set_units(45, "min")`), a `numeric` in seconds, or a string
#'   parseable by `lubridate::duration()` (e.g. `"45 mins"`).
#' @param tolerance Half-width of the acceptance window around
#'   `interval`, in the same units as `interval`. Defaults to
#'   `interval / 10`.
#' @param criterion Tie-breaking rule when multiple chains of equal
#'   length exist: `"closest"` (minimum total deviation from the
#'   target interval, default) or `"first"` (earliest-starting chain).
#' @param remove Logical. If `TRUE`, return only the retained fixes.
#'   Default `FALSE` — the full object is returned with the
#'   `thin_selected` flag.
#'
#' @return A `move2` object. With `remove = FALSE`, the input is
#'   returned unchanged except for a new logical column
#'   `thin_selected`. With `remove = TRUE`, only rows with
#'   `thin_selected == TRUE` are returned. For multi-track input the
#'   thinning is performed independently per track.
#'
#' @examples
#' \donttest{
#' library(move2)
#' fishers <- mt_read(mt_example())
#' fishers <- fishers[!sf::st_is_empty(fishers), ]
#' leroy   <- fishers[mt_track_id(fishers) == "M4", ][seq_len(500), ]
#'
#' out <- mt_thin_time(leroy, interval = "45 min", tolerance = "5 min")
#' table(out$thin_selected)
#'
#' kept <- out[out$thin_selected, ]
#' summary(as.numeric(mt_time_lags(kept, units = "min")))
#' }
#'
#' @seealso [move2::mt_filter_per_interval()] for bucket-aligned
#'   thinning (different semantics), [move2::mt_interpolate()] for the
#'   inverse problem of increasing resolution.
#'
#' @export
mt_thin_time <- function(x,
                         interval,
                         tolerance = NULL,
                         criterion = c("closest", "first"),
                         remove    = FALSE) {
  stopifnot(inherits(x, "move2"))
  criterion <- match.arg(criterion)

  interval_s  <- .parse_duration_secs(interval,  "interval")
  if (is.null(tolerance)) {
    tolerance_s <- interval_s / 10
  } else {
    tolerance_s <- .parse_duration_secs(tolerance, "tolerance")
  }
  if (interval_s <= 0) {
    rlang::abort("`interval` must be positive.",
                 class = "move2utils_mt_thin_time_bad_interval")
  }
  if (tolerance_s < 0) {
    rlang::abort("`tolerance` must be non-negative.",
                 class = "move2utils_mt_thin_time_bad_tolerance")
  }
  if (tolerance_s > interval_s) {
    rlang::abort("`tolerance` must not exceed `interval`.",
                 class = "move2utils_mt_thin_time_tolerance_exceeds_interval")
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
    ts  <- as.numeric(move2::mt_time(x[idx, ]))
    ord <- order(ts)
    ts_sorted <- ts[ord]

    sel_sorted <- .thin_dp_track(ts_sorted, interval_s, tolerance_s,
                                  criterion)
    ## map back to original row order within this track
    sel_in_track <- logical(length(idx))
    sel_in_track[ord] <- sel_sorted
    selected[idx] <- sel_in_track
  }

  x$thin_selected <- selected
  if (remove) x[selected, ] else x
}

## ------------------------------------------------------------------
## Helpers

.parse_duration_secs <- function(value, name) {
  if (inherits(value, "units")) {
    ## A `units` object satisfies is.numeric(); intercept it before the
    ## numeric branch, where `as.numeric()` would strip the unit without
    ## converting (e.g. set_units(45, "min") -> 45 "seconds", not 2700).
    converted <- tryCatch(
      units::set_units(value, "s"),
      error = function(e) {
        rlang::abort(
          sprintf("`%s` is a units object that cannot be converted to seconds (a time unit is required).",
                  name),
          class  = "move2utils_mt_thin_time_bad_units",
          parent = e)
      })
    as.numeric(converted)
  } else if (inherits(value, "difftime")) {
    as.numeric(value, units = "secs")
  } else if (is.numeric(value)) {
    as.numeric(value)
  } else if (is.character(value) && length(value) == 1) {
    if (!requireNamespace("lubridate", quietly = TRUE)) {
      rlang::abort(
        sprintf("String `%s` requires the `lubridate` package.", name),
        class = "move2utils_mt_thin_time_missing_lubridate")
    }
    as.numeric(lubridate::duration(value), units = "secs")
  } else {
    rlang::abort(
      sprintf("`%s` must be a difftime, a numeric in seconds, or a duration string.",
              name),
      class = "move2utils_mt_thin_time_bad_duration_type")
  }
}

## Run-split on gaps > (d + tol), then DP within each run.
.thin_dp_track <- function(ts, d, tol, criterion) {
  n <- length(ts)
  if (n == 0L) return(logical(0))
  if (n == 1L) return(TRUE)

  gaps <- diff(ts)
  run_id <- c(0L, cumsum(gaps > d + tol))
  selected <- logical(n)
  for (r in unique(run_id)) {
    idx <- which(run_id == r)
    selected[idx] <- .thin_dp_run(ts[idx], d, tol, criterion)
  }
  selected
}

## Core DP on a single run of sorted timestamps.
.thin_dp_run <- function(ts, d, tol, criterion) {
  m <- length(ts)
  if (m == 0L) return(logical(0))
  if (m == 1L) return(TRUE)

  lo <- d - tol
  hi <- d + tol

  f    <- rep(1L, m)
  dev  <- rep(0,  m)
  pred <- rep(NA_integer_, m)

  for (i in seq_len(m)[-1]) {
    best_f  <- 1L
    best_j  <- NA_integer_
    best_dv <- Inf
    ti <- ts[i]
    j <- i - 1L
    while (j >= 1L) {
      gap <- ti - ts[j]
      if (gap > hi) break
      if (gap >= lo) {
        fj  <- f[j] + 1L
        dvj <- dev[j] + abs(gap - d)
        take <- FALSE
        if (fj > best_f) {
          take <- TRUE
        } else if (fj == best_f && criterion == "closest" &&
                   dvj < best_dv) {
          take <- TRUE
        }
        if (take) {
          best_f  <- fj
          best_j  <- j
          best_dv <- dvj
        }
      }
      j <- j - 1L
    }
    f[i]    <- best_f
    dev[i]  <- if (is.na(best_j)) 0 else best_dv
    pred[i] <- best_j
  }

  ## pick the end node of the best chain
  max_f <- max(f)
  candidates <- which(f == max_f)
  if (criterion == "closest") {
    end_i <- candidates[which.min(dev[candidates])]
  } else {
    end_i <- candidates[1]              # "first" = earliest ending
  }

  ## backtrace
  selected <- logical(m)
  k <- end_i
  while (!is.na(k)) {
    selected[k] <- TRUE
    k <- pred[k]
  }
  selected
}

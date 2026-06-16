# Corridor detection for move2 objects.
#
# Adapted from the legacy `move::corridor()` function. The deprecated
# `sp`/`rgeos` dependency chain is dropped in favour of `sf` spatial
# indexing (GEOS R-tree). Segment geometry, distance, speed and
# azimuth come from the track-aware `move2` primitives, so the
# function handles multi-track input natively and operates in the
# caller's CRS (longlat is reprojected to a local AEQD only for the
# buffer-and-index step).

#' Identify corridor segments in a move2 trajectory
#'
#' Detects corridor behaviour by finding segments where the animal moves
#' quickly in a spatially consistent direction. Corridors are characterised
#' by high segment speed combined with low circular variance of travel
#' direction among spatially neighbouring segments.
#'
#' The algorithm:
#' \enumerate{
#'   \item Build per-segment linestrings with [move2::mt_segments()] and
#'         take their midpoints with [sf::st_line_interpolate()]. Both
#'         are track-aware; cross-track midpoints are never produced.
#'   \item For each segment, find all other midpoints within its own
#'         half-length search radius via an `sf` spatial index (R-tree).
#'         The search uses the input CRS if projected, or a local
#'         azimuthal-equidistant projection ([move2::mt_aeqd_crs()]) if
#'         the input is in longitude/latitude.
#'   \item Compute a pseudo-azimuth (`(2 * azimuth) mod 2*pi`, working
#'         in radians) so parallel opposite directions collapse onto
#'         the same value: animals walking north and south along the
#'         same corridor are treated as directionally consistent.
#'   \item Compute the circular variance of pseudo-azimuths within each
#'         neighbourhood.
#'   \item Flag a segment as corridor if its speed is at least
#'         `speed_threshold`, its circular variance is at most
#'         `circvar_threshold`, and the corridor-qualifying neighbours
#'         within its search radius outnumber the non-qualifying ones
#'         (at least `min_segments` qualifying neighbours required).
#' }
#'
#' Thresholds are user-supplied. When either is left as `NULL` the
#' function falls back to a within-object quantile (0.75 of segment
#' speeds; 0.25 of valid circular variances) and emits a warning
#' naming the resolved numeric. For comparable corridor maps across
#' individuals or populations, supply explicit thresholds: e.g.
#' `speed_threshold = quantile(as.numeric(mt_speed(all_tracks,
#' units = "m/s")), 0.75, na.rm = TRUE)` computed once on the pooled
#' cohort and reused per individual.
#'
#' This is a port of the `move::corridor()` concept (LaPoint et al. 2013)
#' to the `move2` / `sf` stack, using modern spatial indexing.
#'
#' @param x A `move2` object. Multi-track input is handled track-by-track
#'   via [move2::mt_segments()]; the input may be in any CRS (geographic
#'   or projected). Rows with empty geometries are dropped from the
#'   computation with an informational notice and re-padded with `NA`
#'   in the output.
#' @param speed_threshold Numeric. Per-segment speed (metres / second) at
#'   or above which a segment is "fast." If `NULL` (default), the 0.75
#'   quantile of within-object segment speeds is used and a warning
#'   names the resolved value.
#' @param circvar_threshold Numeric in `[0, 1]`. Circular variance of
#'   pseudo-azimuths within a segment's neighbourhood at or below which
#'   the segment is "directionally consistent." If `NULL` (default),
#'   the 0.25 quantile of valid within-object circular variances is
#'   used and a warning names the resolved value.
#' @param min_segments Integer. Minimum number of qualifying neighbours
#'   required for a segment to be classified as a corridor (default 2).
#' @param verbose Logical. If `TRUE`, prints a one-line summary of how
#'   many segments were flagged. Default `FALSE`. Independent of the
#'   default-threshold, irregular-sampling, and empty-geometry notices,
#'   which are shown regardless.
#'
#' @return The input `move2` object with five added columns, each of
#'   length `nrow(x)`:
#'   \describe{
#'     \item{`corridor`}{factor with levels `"corridor"` and
#'       `"not corridor"`.}
#'     \item{`corridor_speed`}{per-segment speed (m/s).}
#'     \item{`corridor_azimuth`}{per-segment travel azimuth (degrees).}
#'     \item{`corridor_circvar`}{circular variance of pseudo-azimuths
#'       within the neighbourhood.}
#'     \item{`corridor_n_neighbours`}{number of segments found within
#'       the search radius.}
#'   }
#'   Track-final rows, rows with empty geometries, and segments whose
#'   neighbourhood circular variance cannot be computed all carry `NA`
#'   in the per-segment columns; their `corridor` factor level is
#'   `"not corridor"`.
#'
#' @references
#' LaPoint S, Gallery P, Wikelski M, Kays R (2013). Animal behavior,
#' cost-based corridor models, and real corridors. *Landscape Ecology*,
#' 28, 1615–1630. \doi{10.1007/s10980-013-9910-0}
#'
#' @examples
#' \dontrun{
#' library(move2)
#' fishers <- mt_read(mt_example())
#' out <- mt_corridor(fishers)
#' table(out$corridor)
#'
#' ## fixed thresholds for cross-individual comparability:
#' all_speeds <- as.numeric(mt_speed(fishers, units = "m/s"))
#' out2 <- mt_corridor(
#'   fishers,
#'   speed_threshold   = stats::quantile(all_speeds, 0.75, na.rm = TRUE),
#'   circvar_threshold = 0.2
#' )
#' }
#'
#' @export
mt_corridor <- function(x,
                        speed_threshold   = NULL,
                        circvar_threshold = NULL,
                        min_segments      = 2L,
                        verbose           = FALSE) {

  if (!inherits(x, "move2")) {
    rlang::abort(
      "`x` must be a `move2` object.",
      class = "move2utils_input_not_move2"
    )
  }

  n_orig    <- nrow(x)
  is_empty  <- sf::st_is_empty(x)
  n_empty   <- sum(is_empty)
  if (n_empty > 0L) {
    rlang::inform(
      sprintf(
        "Dropping %d row(s) with empty geometries from corridor detection.",
        n_empty
      ),
      class = "move2utils_mt_corridor_empty_dropped"
    )
  }
  keep_idx <- which(!is_empty)
  if (length(keep_idx) < 2L) {
    rlang::abort(
      "Need at least 2 non-empty rows for corridor detection.",
      class = "move2utils_mt_corridor_too_few_rows"
    )
  }
  x_kept <- x[keep_idx, ]

  ## --- track-aware segment geometry -----------------------------------
  seg_geom  <- move2::mt_segments(x_kept)
  geom_type <- as.character(sf::st_geometry_type(seg_geom))
  is_seg    <- geom_type == "LINESTRING"
  n_seg     <- sum(is_seg)
  if (n_seg < 2L) {
    rlang::abort(
      "Need at least 2 segments across all tracks for corridor detection.",
      class = "move2utils_mt_corridor_too_few_segments"
    )
  }

  seg_length_m <- as.numeric(move2::mt_distance(x_kept, units = "m"))[is_seg]
  speed_seg    <- as.numeric(move2::mt_speed(x_kept,    units = "m/s"))[is_seg]
  seg_radius_m <- seg_length_m / 2
  ## Azimuth is computed in the work CRS (Euclidean atan2 on segment
  ## endpoints) further below, after `seg_lines` has been transformed
  ## to a metric CRS.  This keeps the function projection-agnostic
  ## without the longlat detour that move2::mt_azimuth() would force.

  ## --- sampling-interval sanity warning -------------------------------
  ## Half-segment-length search radius implicitly assumes roughly
  ## uniform sampling. On wildly irregular tracks, segment lengths
  ## (and so search radii) span orders of magnitude and corridor flags
  ## on the long-tail segments become hard to interpret.
  ## See HEURISTICS.md (Group 2) for the IQR/median ratio rationale.
  time_lags_num <- suppressWarnings(
    as.numeric(move2::mt_time_lags(x_kept, units = "s"))
  )
  tl_finite <- time_lags_num[is.finite(time_lags_num) & time_lags_num > 0]
  if (length(tl_finite) >= 2L) {
    tl_med <- stats::median(tl_finite)
    tl_iqr <- stats::IQR(tl_finite)
    if (tl_med > 0 && tl_iqr / tl_med > 1) {
      rlang::warn(
        sprintf(
          paste0(
            "Sampling interval is highly irregular ",
            "(IQR / median = %.2g). The half-segment-length search ",
            "radius assumes roughly uniform sampling; corridor flags ",
            "on very short or very long segments may be misleading."
          ),
          tl_iqr / tl_med
        ),
        class = "move2utils_mt_corridor_irregular_sampling"
      )
    }
  }

  ## --- midpoints + neighbour search + Euclidean azimuth ---------------
  ## st_buffer + st_intersects need a metric CRS, and st_line_interpolate
  ## with `normalized = TRUE` is only meaningful on planar input.  Transform
  ## segment lines to a local AEQD via mt_aeqd_crs() when the caller's CRS
  ## is longlat; otherwise trust the user's projected CRS.  Azimuth is then
  ## the planar bearing atan2(dx, dy) computed directly on the transformed
  ## segment endpoints — Cartesian on a projected CRS, AEQD-planar on
  ## longlat input.  AEQD distortion is sub-degree for typical
  ## movement-ecology segment lengths.
  seg_lines <- seg_geom[is_seg]
  if (isTRUE(sf::st_is_longlat(x_kept))) {
    work_crs  <- move2::mt_aeqd_crs(x_kept)
    seg_lines <- sf::st_transform(seg_lines, crs = work_crs)
  }
  mid_proj   <- sf::st_line_interpolate(seg_lines, dist = 0.5,
                                        normalized = TRUE)
  buffered   <- sf::st_buffer(mid_proj, dist = seg_radius_m)
  neighbours <- sf::st_intersects(buffered, mid_proj)
  n_neighbours <- lengths(neighbours)

  ## Endpoint coordinates: mt_segments() emits 2-point LINESTRINGs,
  ## so st_coordinates() returns row pairs (start, end) per segment.
  seg_cc      <- sf::st_coordinates(seg_lines)
  start_xy    <- seg_cc[seq(1L, by = 2L, length.out = n_seg), c("X", "Y")]
  end_xy      <- seg_cc[seq(2L, by = 2L, length.out = n_seg), c("X", "Y")]
  azimuth_rad <- atan2(end_xy[, "X"] - start_xy[, "X"],
                       end_xy[, "Y"] - start_xy[, "Y"])
  ## Pseudo-azimuth: doubling collapses opposite-direction segments
  ## onto the same value, so two animals walking the same corridor in
  ## opposite directions are treated as directionally consistent.
  ## Formula in radians: pseudo = (2 * az) mod 2*pi.
  pseudo_az_rad <- (2 * azimuth_rad) %% (2 * pi)

  ## --- circular variance per neighbourhood ----------------------------
  circ_var <- .circ_var_per_neighbourhood(pseudo_az_rad, neighbours)

  ## --- threshold resolution -------------------------------------------
  if (is.null(speed_threshold)) {
    speed_threshold <- stats::quantile(speed_seg, probs = 0.75, na.rm = TRUE)
    rlang::warn(
      sprintf(
        paste0(
          "`speed_threshold` not supplied; using within-object 0.75 ",
          "quantile of segment speed = %.4g m/s. For cross-individual ",
          "or cross-population comparability, supply an explicit value."
        ),
        as.numeric(speed_threshold)
      ),
      class = "move2utils_mt_corridor_default_speed_threshold"
    )
  }
  if (is.null(circvar_threshold)) {
    cv_valid <- circ_var[!is.na(circ_var)]
    circvar_threshold <- if (length(cv_valid) >= 2L) {
      stats::quantile(cv_valid, probs = 0.25, na.rm = TRUE)
    } else {
      NA_real_
    }
    rlang::warn(
      sprintf(
        paste0(
          "`circvar_threshold` not supplied; using within-object 0.25 ",
          "quantile of valid circular variance = %.4g. For ",
          "cross-individual or cross-population comparability, supply ",
          "an explicit value."
        ),
        as.numeric(circvar_threshold)
      ),
      class = "move2utils_mt_corridor_default_circvar_threshold"
    )
  }

  is_fast       <- !is.na(speed_seg) & speed_seg >= as.numeric(speed_threshold)
  is_consistent <- !is.na(circ_var) &
                   circ_var <= as.numeric(circvar_threshold)
  is_candidate  <- is_fast & is_consistent

  ## --- candidate-majority test per neighbourhood ----------------------
  is_corridor <- rep(FALSE, n_seg)
  if (any(is_candidate)) {
    owner      <- rep.int(seq_len(n_seg), n_neighbours)
    nb_flat    <- unlist(neighbours, use.names = FALSE)
    own_f      <- factor(owner, levels = seq_len(n_seg))
    cand_in_nb <- is_candidate[nb_flat]
    n_cand     <- as.numeric(tapply(cand_in_nb, own_f, sum))
    is_corridor <- n_neighbours >= min_segments &
                   n_cand > (n_neighbours - n_cand)
  }

  ## --- assemble length-n_kept output, then NA-pad to length-n_orig ----
  n_kept <- nrow(x_kept)
  seg_speed_kept   <- rep(NA_real_,    n_kept)
  seg_azimuth_kept <- rep(NA_real_,    n_kept)
  seg_circvar_kept <- rep(NA_real_,    n_kept)
  seg_nb_kept      <- rep(NA_integer_, n_kept)
  corridor_label_kept <- rep("not corridor", n_kept)

  seg_speed_kept[is_seg]   <- speed_seg
  ## store azimuth in degrees (matching the @return docstring).
  seg_azimuth_kept[is_seg] <- azimuth_rad * 180 / pi
  seg_circvar_kept[is_seg] <- circ_var
  seg_nb_kept[is_seg]      <- as.integer(n_neighbours)
  corridor_label_kept[which(is_seg)[is_corridor]] <- "corridor"

  corridor_label   <- rep("not corridor", n_orig)
  corridor_speed   <- rep(NA_real_,    n_orig)
  corridor_azimuth <- rep(NA_real_,    n_orig)
  corridor_circvar <- rep(NA_real_,    n_orig)
  corridor_nb      <- rep(NA_integer_, n_orig)
  corridor_label[keep_idx]   <- corridor_label_kept
  corridor_speed[keep_idx]   <- seg_speed_kept
  corridor_azimuth[keep_idx] <- seg_azimuth_kept
  corridor_circvar[keep_idx] <- seg_circvar_kept
  corridor_nb[keep_idx]      <- seg_nb_kept

  x$corridor              <- factor(corridor_label,
                                    levels = c("corridor", "not corridor"))
  x$corridor_speed        <- corridor_speed
  x$corridor_azimuth      <- corridor_azimuth
  x$corridor_circvar      <- corridor_circvar
  x$corridor_n_neighbours <- corridor_nb

  if (isTRUE(verbose)) {
    n_corr <- sum(is_corridor)
    message(sprintf(
      "Found %d corridor segments (%.1f%% of %d segments).",
      n_corr, 100 * n_corr / n_seg, n_seg
    ))
  }
  x
}

## Internal: vectorised circular variance per neighbourhood.
##
## For pseudo-azimuths in radians, circular variance is 1 - |R| where
## R = (cos_mean, sin_mean) is the mean resultant vector. We flatten
## the neighbour list once and aggregate sin/cos sums per owner in a
## single pass instead of calling circular::var.circular() per segment.
## Returns numeric of length(neighbour_list); NA for neighbourhoods of
## size < 2.
.circ_var_per_neighbourhood <- function(angles_rad, neighbour_list) {
  n      <- length(neighbour_list)
  nl     <- lengths(neighbour_list)
  cos_a  <- cos(angles_rad)
  sin_a  <- sin(angles_rad)
  owner  <- rep.int(seq_len(n), nl)
  nb_flat <- unlist(neighbour_list, use.names = FALSE)
  own_f  <- factor(owner, levels = seq_len(n))
  sum_cos <- as.numeric(tapply(cos_a[nb_flat], own_f, sum))
  sum_sin <- as.numeric(tapply(sin_a[nb_flat], own_f, sum))
  R <- sqrt((sum_cos / nl)^2 + (sum_sin / nl)^2)
  cv <- 1 - R
  cv[nl < 2L] <- NA_real_
  cv
}

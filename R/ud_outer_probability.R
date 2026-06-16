#' Outer probability at a set of locations
#'
#' Given a utilisation distribution and a set of query points, return
#' for each point the **cumulative-volume quantile** at that location:
#' the fraction of the UD contained in cells of equal or higher density.
#' If a query point falls on a cell whose value is 0.5 in the volume
#' UD, the point is on the 50 %-isopleth — together with more-central
#' cells, those cells hold half of the total distribution.
#'
#' Equivalent to asking "how peripheral is this location within the
#' animal's space use?". A value near 0 means the location is at a
#' density peak; a value near 1 means it is far out in the tail.
#'
#' This is the move2 analogue of `move::outerProbability()`, and is a
#' thin wrapper around [ud_volume()] + [terra::extract()].
#'
#' @param x An `sf` or `move2` object containing the query points, or
#'   a two-column numeric matrix/data-frame of `x,y` coordinates. `sf`
#'   and `move2` query points are reprojected to the CRS of `ud`
#'   automatically; a bare coordinate matrix carries no CRS and must
#'   therefore already be in the same CRS as `ud`.
#' @param ud A single-layer `terra::SpatRaster` holding a
#'   probability-density UD (not a volume UD — the function computes
#'   the volume transform internally).
#' @param ... Additional arguments passed to [terra::extract()].
#'
#' @return A numeric vector with one value per row of `x`: the volume
#'   quantile at each query location. `NA` where the location falls
#'   outside the raster or on an `NA` cell.
#'
#' @examples
#' \dontrun{
#' library(move2)
#' library(sf)
#' fishers <- mt_read(mt_example())
#' fishers <- fishers[!st_is_empty(fishers), ]
#' leroy <- fishers[mt_track_id(fishers) == "M4", ][seq_len(200), ]
#' leroy <- st_transform(leroy, mt_aeqd_crs(leroy))
#'
#' ud <- mt_dbbmm_ud(leroy, location_error = 20, raster = 100)
#'
#' ## outer probability at every fix along the track
#' op <- ud_outer_probability(leroy, ud)
#' summary(op)
#' }
#'
#' @seealso [ud_volume()] for the underlying transform;
#'   [mt_dbbmm_ud()] / [mt_dbgb_ud()] for producing UDs.
#'
#' @references
#' Fieberg, J., & Kochanny, C. O. (2005). Quantifying home-range
#' overlap: the importance of the utilization distribution. *Journal
#' of Wildlife Management*, 69(4), 1346-1359.
#' \doi{10.2193/0022-541X(2005)69[1346:QHOTIO]2.0.CO;2}
#'
#' @export
ud_outer_probability <- function(x, ud, ...) {
  if (!inherits(ud, "SpatRaster")) {
    rlang::abort("`ud` must be a terra::SpatRaster.",
                 class = "move2utils_ud_outer_probability_ud_not_spatraster")
  }
  if (terra::nlyr(ud) != 1L) {
    rlang::abort("`ud` must be a single-layer SpatRaster.",
                 class = "move2utils_ud_outer_probability_ud_multi_layer")
  }

  if (inherits(x, c("sf", "sfc", "move2"))) {
    ## Align query points to the UD's CRS.  Since UDs are returned in the
    ## caller's CRS (which may differ from where the query points live),
    ## extracting raw coordinates against the raster would silently return
    ## NA / wrong cells on a CRS mismatch.  Reprojecting points is exact.
    ud_crs <- sf::st_crs(terra::crs(ud))
    x_crs  <- sf::st_crs(x)
    if (!is.na(ud_crs) && !is.na(x_crs) && ud_crs != x_crs) {
      x <- sf::st_transform(x, ud_crs)
    }
    xy <- sf::st_coordinates(x)[, c(1, 2), drop = FALSE]
  } else if (is.matrix(x) || is.data.frame(x)) {
    if (ncol(x) < 2) {
      rlang::abort("Numeric `x` must have at least two columns for x, y.",
                   class = "move2utils_ud_outer_probability_bad_matrix_cols")
    }
    xy <- as.matrix(x)[, c(1, 2), drop = FALSE]
  } else {
    rlang::abort(paste0(
      "`x` must be an sf / sfc / move2 object, or a 2-column ",
      "coordinate matrix."),
      class = "move2utils_ud_outer_probability_bad_x_type")
  }

  vud <- ud_volume(ud)
  as.numeric(terra::extract(vud, xy, ...)[[1]])
}

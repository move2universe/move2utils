#' Convert a utilisation distribution to cumulative-volume space
#'
#' Transforms a probability-density utilisation distribution (UD) into
#' its cumulative-volume representation: each cell's value becomes the
#' fraction of the distribution contained in cells of equal or higher
#' density. A cell value of *v* in the output means "this cell is
#' inside the *v*-volume contour, together with all denser cells".
#'
#' This is the move2 analogue of `move::getVolumeUD()`. It is a
#' reshape of the raster only — no re-estimation — so it is cheap
#' and can be applied liberally to already-computed UDs.
#'
#' @param x A single-layer or multi-layer `terra::SpatRaster`
#'   representing one or more UDs. Cells should be non-negative; the
#'   UD need not sum exactly to 1 (it is renormalised per layer
#'   internally).
#'
#' @return A `terra::SpatRaster` with the same geometry as `x`, cell
#'   values replaced by their cumulative-volume quantile, and layer
#'   names preserved. NA cells stay NA.
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
#' ud  <- mt_dbbmm_ud(leroy, location_error = 20, raster = 100)
#' vud <- ud_volume(ud)
#'
#' ## draw the 50% and 95% contours of the volume UD
#' terra::plot(vud)
#' terra::contour(vud, levels = c(0.5, 0.95), add = TRUE)
#'
#' ## solid 50% core-area mask
#' core50 <- vud <= 0.5
#' }
#'
#' @seealso [mt_dbbmm_ud()], [mt_dbgb_ud()] for computing UDs;
#'   [terra::contour()] for isopleth extraction from a volume UD.
#'
#' @references
#' Seaman, D. E., & Powell, R. A. (1996). An evaluation of the accuracy
#' of kernel density estimators for home range analysis. *Ecology*,
#' 77(7), 2075-2085. \doi{10.2307/2265701}
#'
#' @export
ud_volume <- function(x) {
  if (!inherits(x, "SpatRaster")) {
    rlang::abort("`x` must be a terra::SpatRaster.",
                 class = "move2utils_ud_volume_not_spatraster")
  }

  out <- x
  for (i in seq_len(terra::nlyr(x))) {
    v <- terra::values(x[[i]])[, 1]
    vtot <- sum(v, na.rm = TRUE)
    if (!is.finite(vtot) || vtot <= 0) {
      next
    }
    v_norm <- v / vtot
    ok <- !is.na(v_norm)
    ord <- order(v_norm, decreasing = TRUE)
    cum <- numeric(length(v_norm))
    cum[ord] <- cumsum(v_norm[ord])
    cum[!ok] <- NA_real_
    terra::values(out[[i]]) <- cum
  }
  out
}

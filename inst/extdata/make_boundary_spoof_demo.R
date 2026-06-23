## ==========================================================================
## Boundary-spoof block-expansion fixture for move2utils.
##
## WHY THIS EXISTS: the shipped synthetic cohort (CPF_A--F) does NOT exercise
## graph-based block expansion. CPF_D's spoof block sits *mid-track*, so cutting
## the kept-fix graph at its impossible boundary transitions severs the
## legitimate trajectory into two comparable halves and the dominance gate
## correctly declines; the block is recovered by the per-fix/consensus layer
## instead. To exercise (and regression-test) block expansion you need a
## *sustained, boundary-anchored* block, which this generator builds.
##
## Construction: the clean reference track CPF_B (n = 3537) with its final 150
## fixes overwritten by a coherent spoof block pinned ~120 km off-route (a
## physiologically impossible boundary transition) with tight within-block
## jitter, so the block interior looks locally normal. Deterministic (seeded).
##
## RETURNS a list:
##   $track : move2 object (EPSG:4326), the boundary-spoof track
##   $truth : integer row indices of the injected spoof block
##   $block : alias of $truth
##
## USAGE:
##   source(system.file("extdata", "make_boundary_spoof_demo.R", package = "move2utils"))
##   d  <- make_boundary_spoof_demo()
##   o  <- mt_clean_track(d$track)                          # auto path
##   sum(o$is_outlier[d$truth])                             # 150 (block-expansion fires)
##   sum(!is.na(o$block_id))                                # ~146 (error_class == "block")
##   oc <- mt_clean_track(d$track, mass = 1, mode = "flying")
##   sum(oc$is_outlier[d$truth])                            # 150 -- a supplied cap recovers
##                                                          #   the block just as well (it is
##                                                          #   the connectivity ceiling block
##                                                          #   expansion uses)
##   sum(!is.na(oc$block_id))                               # ~148
##
## Recovery is consistent across the auto and `mass`/`mode` paths by design:
## supplying a physiological cap never reduces block recovery.  Regression-
## tested in tests/testthat/test-mt_clean_track.R.
## ==========================================================================

make_boundary_spoof_demo <- function(block = 150L, seed = 42L,
                                     offset_deg = 1.5, jitter_sd = 5e-4) {
  stopifnot(requireNamespace("move2", quietly = TRUE),
            requireNamespace("sf", quietly = TRUE))
  m <- utils::read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                          package = "move2utils")),
                       stringsAsFactors = FALSE)
  m$timestamp <- as.POSIXct(m$timestamp, tz = "UTC")
  b <- m[m$individual.local.identifier == "CPF_B", ]
  b <- b[order(b$timestamp), ]
  b$individual.local.identifier <- "SPOOF_BOUNDARY"

  n <- nrow(b)
  idx <- (n - block + 1L):n              # block at the END (boundary, not mid-track)
  set.seed(seed)
  b$location.long[idx] <- b$location.long[1] + offset_deg + stats::rnorm(block, 0, jitter_sd)
  b$location.lat[idx]  <- b$location.lat[1]  + offset_deg + stats::rnorm(block, 0, jitter_sd)

  trk <- move2::mt_as_move2(b, coords = c("location.long", "location.lat"),
                            time_column = "timestamp",
                            track_id_column = "individual.local.identifier",
                            crs = 4326)
  trk <- trk[!sf::st_is_empty(trk), ]
  trk <- trk[order(move2::mt_time(trk)), ]

  ## CPF_B has no empty geometries, so the spoof block remains the final
  ## `block` rows after filtering; compute truth indices on the final object.
  truth <- (nrow(trk) - block + 1L):nrow(trk)
  list(track = trk, truth = truth, block = truth)
}

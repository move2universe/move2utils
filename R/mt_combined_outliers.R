#' Combined outlier detection using multiple methods
#'
#' Runs the simultaneous method (with gap and entropy thresholds) and
#' the sequential scan, then combines their votes.  A location flagged
#' by at least \code{min_votes} methods is declared an outlier.
#'
#' @details
#' \strong{Relationship to the four-primitive cascade.}  This function
#' votes across three \emph{strategies} on a single detector (the
#' joint-probability surface from \code{\link{mt_flag_outliers}}):
#' simultaneous-with-gap, simultaneous-with-entropy, and the
#' sequential scan from \code{\link{mt_sequential_outliers}}.  The
#' unified cleaner \code{\link{mt_clean_track}} is a different
#' construction: it votes across four \emph{detectors} (bridge,
#' detour, probability, speed-cap), each looking at a different
#' grain of the data.  The two are not redundant -- the cascade
#' cannot reach inside a single detector to compare simultaneous vs
#' sequential evaluation.  Reach for \code{mt_combined_outliers}
#' when you want to inspect the per-strategy agreement on the
#' probability surface specifically, or when calibrating thresholds
#' on that surface and wanting cross-strategy validation.  For
#' routine cleaning, \code{mt_clean_track} is the recommended
#' entry-point.
#'
#' @param x A \code{move2} object.
#' @param reference Optional \code{move2} object for reference
#'   distributions.
#' @param min_votes Minimum number of methods (out of 3) that must flag
#'   a location for it to be declared an outlier.  Default: 2 (majority
#'   vote).
#' @param scan Scanning strategy for the sequential method:
#'   \code{"forward-backward"}, \code{"greedy"}, or \code{"random"}.
#' @param time_normalize Logical; if TRUE (default), use speed and
#'   angular velocity.
#' @param plot Logical; if TRUE, plot the results.
#' @param ... Additional arguments passed to \code{mt_flag_outliers} and
#'   \code{mt_sequential_outliers}.
#'
#' @return The input \code{move2} object with added columns:
#'   \code{is_outlier}, \code{vote_count} (0--3), and the individual
#'   method flags \code{flag_gap}, \code{flag_entropy}, \code{flag_seq}.
#'
#' @references
#' Safi, K. (2026). Self-thresholding hierarchical
#' outlier-detection for animal movement tracks. Companion paper to
#' the \pkg{move2utils} R package. bioRxiv preprint, submitted to
#' Methods in Ecology and Evolution. \doi{10.64898/2026.07.11.737894}
#'
#' @seealso \code{\link{mt_clean_track}} (recommended unified
#'   cleaner; votes across four \emph{detectors} rather than three
#'   strategies on one detector); \code{\link{mt_flag_outliers}}
#'   (the probability primitive this function votes strategies on);
#'   \code{\link{mt_sequential_outliers}} (one of the three voted
#'   strategies); \code{\link{mt_persistence_score}} (multi-scale
#'   persistence annotation that can be applied to this function's
#'   output for additional confidence quantification).
#'
#' @examples
#' \dontrun{
#' ## Majority-vote flagging across three detection strategies:
#' res <- mt_combined_outliers(track, min_votes = 2,
#'                              scan = "forward-backward")
#' table(res$vote_count)
#' }
#'
#' @export
mt_combined_outliers <- function(x, reference = NULL,
                                 min_votes = 2,
                                 scan = "forward-backward",
                                 time_normalize = TRUE,
                                 plot = FALSE, ...) {
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  .reject_empty_geometry(x, "mt_combined_outliers()")

  message("=== Combined outlier detection (3-method voting) ===")

  ## ---- method 1: simultaneous with gap threshold ----
  message("\n--- Method 1: simultaneous (gap threshold) ---")
  r_gap <- suppressMessages(
    mt_flag_outliers(x, plot = FALSE, threshold_type = "gap",
                     time_normalize = time_normalize, ...)
  )
  flag_gap <- r_gap$is_outlier
  message(sprintf("  Gap: %d flagged", sum(flag_gap, na.rm = TRUE)))

  ## ---- method 2: simultaneous with entropy threshold ----
  message("--- Method 2: simultaneous (entropy threshold) ---")
  r_ent <- suppressMessages(
    mt_flag_outliers(x, plot = FALSE, threshold_type = "entropy",
                     time_normalize = time_normalize, ...)
  )
  flag_ent <- r_ent$is_outlier
  message(sprintf("  Entropy: %d flagged", sum(flag_ent, na.rm = TRUE)))

  ## ---- method 3: sequential scan ----
  message("--- Method 3: sequential scan ---")
  r_seq <- suppressMessages(
    mt_sequential_outliers(x, reference = reference, scan = scan,
                            time_normalize = time_normalize, ...)
  )
  flag_seq <- r_seq$is_outlier
  message(sprintf("  Sequential: %d flagged", sum(flag_seq, na.rm = TRUE)))

  ## ---- combine votes ----
  votes <- as.integer(flag_gap) + as.integer(flag_ent) + as.integer(flag_seq)
  is_outlier <- votes >= min_votes

  ## ---- attach results ----
  x$flag_gap     <- flag_gap
  x$flag_entropy <- flag_ent
  x$flag_seq     <- flag_seq
  x$vote_count   <- votes
  x$is_outlier   <- is_outlier

  ## carry forward the simultaneous probabilities for diagnostics
  x$joint_prob     <- r_gap$joint_prob
  x$seq_joint_prob <- r_seq$seq_joint_prob

  n_out <- sum(is_outlier, na.rm = TRUE)
  message(sprintf(
    "\n=== Combined: %d outliers (min_votes=%d, %.2f%% of %d locations) ===",
    n_out, min_votes, 100 * n_out / nrow(x), nrow(x)
  ))
  message(sprintf("  By vote count: 0=%d, 1=%d, 2=%d, 3=%d",
    sum(votes == 0, na.rm = TRUE), sum(votes == 1, na.rm = TRUE),
    sum(votes == 2, na.rm = TRUE), sum(votes == 3, na.rm = TRUE)
  ))

  if (plot) {
    ids <- as.character(move2::mt_track_id(x))
    unique_ids <- unique(ids)
    n_tracks <- length(unique_ids)
    coords <- sf::st_coordinates(x)
    vote_col <- c("grey80", "#FDAE6B", "#E6550D", "#8B0000")

    ## Multi-track: one panel per individual so tracks aren't visually
    ## overlaid.  Single-track: keep the original single-panel layout.
    if (n_tracks > 1L) {
      op <- graphics::par(mfrow = grDevices::n2mfrow(n_tracks))
      on.exit(graphics::par(op), add = TRUE)
    }

    for (uid in unique_ids) {
      sel <- which(ids == uid)
      main_label <- if (n_tracks > 1L) {
        sprintf("%s (min votes: %d)", uid, min_votes)
      } else {
        paste("Combined outlier detection (min votes:", min_votes, ")")
      }
      graphics::plot(coords[sel, , drop = FALSE], type = "n", asp = 1,
                     xlab = "Longitude", ylab = "Latitude",
                     main = main_label)
      for (v in 0:3) {
        pts <- sel[which(votes[sel] == v)]
        if (length(pts) > 0L) {
          graphics::points(coords[pts, , drop = FALSE],
                           col = vote_col[v + 1], pch = 19,
                           cex = 0.3 + v * 0.3)
        }
      }
      graphics::legend("topright",
                       legend = paste0(0:3, " votes"),
                       col = vote_col, pch = 19,
                       pt.cex = 0.3 + (0:3) * 0.3,
                       bty = "n")
    }
  }

  x
}

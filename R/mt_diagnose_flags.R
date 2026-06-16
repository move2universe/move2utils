#' Per-detector flag audit for a cleaned move2 object
#'
#' Reads the per-detector flag columns that \code{\link{mt_clean_track}}
#' (or any of the four primitives) attaches to its output, returns four
#' diagnostic tables plus a consensus-mode comparison, and a map
#' showing where each detector fired versus where the consensus rule
#' produced the final \code{is_outlier} decision.  Use it to decide
#' whether the cascade is doing what you expect, and where to tune.
#'
#' @details
#' The four tables answer four questions:
#' \enumerate{
#'   \item \strong{error_class} -- where in the cascade hierarchy do
#'     the flags land?  A flag-heavy
#'     \code{kinematic_confluence} bucket means the per-fix detectors
#'     are agreeing 2-of-3.  A heavy \code{block} bucket means block-
#'     expansion is doing the catching.  A heavy \code{pool} bucket
#'     means cohort \code{pool_by} is supplying flags the per-track
#'     cascade missed.
#'   \item \strong{per-detector fires among flagged} -- who is doing
#'     the work?  If one detector fires on >80\% of flags and another
#'     fires on <20\%, your data plays well to one type of evidence and
#'     poorly to another.
#'   \item \strong{co-fire histogram} -- how many detectors typically
#'     agree?  Most flagged fixes should have 2 or more agreeing.
#'     Many one-detector fires that the consensus rejected means the
#'     loudest voice is being silenced.
#'   \item \strong{near-miss zone} -- per detector, how many fixes
#'     fired but were NOT flagged because no other detector
#'     corroborated?  A detector with thousands of near-misses is the
#'     candidate for either tighter thresholds (so other detectors
#'     catch up) or a looser consensus rule (so its verdict is
#'     accepted more often).
#' }
#'
#' The \strong{consensus comparison} re-applies the consensus rule
#' alone (no re-running the per-fix detectors) under each built-in
#' mode of \code{\link{mt_flag_consensus}} and reports how many
#' fixes each mode would flag.  This is the cheap interactive piece:
#' switching consensus mode is post-hoc, so you see the cost of every
#' choice without paying the cascade's compute cost.
#'
#' The \strong{map} has five panels.  Top row: four small panels, one
#' per detector, showing where that detector fired.  Within each
#' panel, fires that survived the consensus are filled, fires that
#' did NOT survive (the near-miss zone) are hollow.  Bottom row: one
#' large panel showing the kept trajectory plus the flagged fixes
#' coloured by \code{error_class}.  All panels share one azimuthal-
#' equidistant projection centred on the object's centroid, so
#' across-detector visual comparison is honest.
#'
#' @param x A \code{move2} object with the per-detector flag columns
#'   that \code{\link{mt_clean_track}} or the \code{mt_flag_*}
#'   primitives attach.  Must contain at minimum \code{is_outlier},
#'   plus a subset of \code{flagged_by_bridge}, \code{flagged_by_prob},
#'   \code{flagged_by_detour}, \code{flagged_by_speed}.
#'   \code{error_class} is required for the bottom map panel and is
#'   produced by \code{\link{mt_clean_track}}.
#' @param print_tables Logical.  If \code{TRUE} (default), print the
#'   four diagnostic tables plus the consensus-comparison table to
#'   the console.  Set \code{FALSE} to silence the console output and
#'   work from the returned list.
#'
#' @section Comparing diagnostics across pool_by configurations:
#' The diagnostic is most informative when run on the same data
#' with different \code{pool_by} settings, to evaluate whether
#' cohort-level pooling is doing useful work or pushing the
#' cleaning past outlier-detection into within-animal behavioural
#' anomaly territory.
#'
#' Pattern: re-run \code{\link{mt_clean_track}} with a narrower
#' \code{pool_by} (e.g. \code{pool_by = "individual_id"} instead of
#' \code{pool_by = c("study_id", "individual_id")}) and compare the
#' two diagnostic outputs.  Three readings to look for:
#'
#' \itemize{
#'   \item \strong{Same total flag count and same per-detector fire
#'     pattern}: the cohort distribution and the per-individual
#'     distribution agree; the pool_by choice is making no
#'     difference and either is fine.
#'   \item \strong{Narrower pool_by flags more on some individuals}:
#'     those animals have movement patterns narrower than the
#'     cohort.  The extra flags may be real outliers the cohort
#'     distribution missed (in which case narrower is better) OR
#'     normal-for-this-animal behaviour that fell outside its own
#'     narrow distribution (in which case broader is better and the
#'     extra flags are false positives).
#'   \item \strong{Broader pool_by flags more}: the cohort has a
#'     heavier tail than the individuals (unusual; usually means a
#'     mix of distinct movement modes across animals).  Cohort
#'     pooling is then aggregating the wrong thing; consider
#'     stratifying by movement type before pooling.
#' }
#'
#' Visual check is decisive when narrower pooling flags more.  Map
#' the newly-flagged fixes; check their step speeds against the
#' physiological cap.  Speeds well below physiology + spatial
#' position inside the kept cluster = within-animal behavioural
#' variability, not outliers.  Speeds near or above physiology +
#' spatial isolation from the home-range bulk = real outliers the
#' broader pooling missed.
#'
#' On datasets where the cohort step-speed distribution is
#' unimodal and the species shows a single dominant movement mode,
#' cohort-level pooling is usually the principled outlier-detection
#' choice.  Narrower pooling crosses into behavioural-anomaly
#' detection, which is a different operation than removing
#' measurement errors.
#'
#' @return Invisibly, a list with components:
#'   \describe{
#'     \item{\code{error_class}}{Frequency table of \code{error_class}
#'       among flagged fixes.}
#'     \item{\code{detector_fires}}{Data frame: for each detector
#'       column present, the count and percentage of flagged fixes
#'       that detector fired on.}
#'     \item{\code{co_fire}}{Two-way table of (n detectors fired) x
#'       (is_outlier) over all events.}
#'     \item{\code{near_miss}}{Data frame: per detector, the number
#'       of fixes where it fired but \code{is_outlier} stayed
#'       \code{FALSE} (consensus did not trip).}
#'     \item{\code{consensus_comparison}}{Data frame: under each
#'       built-in consensus mode (\code{"class_aware"}, \code{"strict"},
#'       \code{"majority"}, \code{"speed_trusted"}, \code{"any"}),
#'       how many fixes would be flagged, expressed in absolute
#'       counts and as a delta from the current \code{is_outlier}
#'       column.}
#'     \item{\code{map}}{A \code{patchwork} object combining the five
#'       map panels.  Use \code{print(result$map)} to draw it.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(move2)
#' x <- mt_read(system.file("extdata/synthetic_tracks.csv.gz",
#'                            package = "move2utils"))
#' x <- x[!sf::st_is_empty(x), ]
#' cleaned <- mt_clean_track(x, plot = FALSE, remove = FALSE)
#' diag <- mt_diagnose_flags(cleaned)
#' print(diag$map)
#' diag$consensus_comparison
#' }
#'
#' @seealso \code{\link{mt_clean_track}} for the orchestrator that
#'   attaches the flag columns; \code{\link{mt_flag_consensus}} for
#'   the post-hoc consensus rule the comparison column re-applies.
#'
#' @importFrom move2 mt_track_id mt_aeqd_crs
#' @importFrom sf st_transform st_coordinates
#' @export
mt_diagnose_flags <- function(x, print_tables = TRUE) {
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object carrying mt_clean_track flag columns.",
                 class = "move2utils_input_not_move2")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    rlang::abort(paste0(
      "Package 'ggplot2' is required for mt_diagnose_flags(). ",
      "Install with install.packages('ggplot2')."),
      class = "move2utils_mt_diagnose_flags_missing_ggplot2")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    rlang::abort(paste0(
      "Package 'patchwork' is required for mt_diagnose_flags(). ",
      "Install with install.packages('patchwork')."),
      class = "move2utils_mt_diagnose_flags_missing_patchwork")
  }
  if (!("is_outlier" %in% names(x))) {
    rlang::abort(paste0(
      "`x` is missing `is_outlier`.  Call mt_clean_track() with ",
      "`remove = FALSE` and pass the result here."),
      class = "move2utils_mt_diagnose_flags_missing_is_outlier")
  }

  fire_candidates <- c("flagged_by_bridge", "flagged_by_prob",
                       "flagged_by_speed", "flagged_by_detour")
  fire_cols <- intersect(fire_candidates, names(x))
  if (!length(fire_cols)) {
    rlang::abort(paste0(
      "`x` has no per-detector flag columns (flagged_by_*).  This ",
      "diagnostic only applies to output of mt_clean_track() or the ",
      "mt_flag_* primitives."),
      class = "move2utils_mt_diagnose_flags_missing_flag_columns")
  }

  ev <- as.data.frame(x)
  is_out <- as.logical(ev$is_outlier)
  is_out[is.na(is_out)] <- FALSE

  ## ---- Table 1: error_class breakdown -------------------------------
  table_error_class <- if ("error_class" %in% names(ev)) {
    as.data.frame(table(error_class = ev$error_class[is_out],
                        useNA = "ifany"),
                  stringsAsFactors = FALSE)
  } else {
    data.frame(error_class = character(), Freq = integer())
  }

  ## ---- Table 2: per-detector fires among flagged ---------------------
  n_flag <- sum(is_out)
  table_detector_fires <- data.frame(
    detector = fire_cols,
    n_fires  = vapply(fire_cols, function(c) sum(ev[[c]] & is_out, na.rm = TRUE),
                       integer(1)),
    stringsAsFactors = FALSE)
  table_detector_fires$pct_of_flagged <- if (n_flag == 0) {
    rep(0, nrow(table_detector_fires))
  } else {
    round(100 * table_detector_fires$n_fires / n_flag, 1)
  }

  ## ---- Table 3: co-fire histogram -----------------------------------
  fired_mat <- vapply(fire_cols, function(c) as.integer(ev[[c]]),
                       integer(nrow(ev)))
  fired_mat[is.na(fired_mat)] <- 0L
  n_fired <- rowSums(fired_mat)
  table_co_fire <- as.data.frame(
    table(n_detectors_fired = n_fired, is_outlier = is_out),
    stringsAsFactors = FALSE)

  ## ---- Table 4: near-miss zone --------------------------------------
  near_miss <- !is_out & (n_fired >= 1)
  table_near_miss <- data.frame(
    detector  = fire_cols,
    n_near_miss = vapply(fire_cols,
                          function(c) sum(ev[[c]] & near_miss, na.rm = TRUE),
                          integer(1)),
    stringsAsFactors = FALSE)

  ## ---- Table 5: consensus comparison --------------------------------
  ## Re-apply each built-in consensus mode (post-hoc, no detector
  ## re-run).  mt_flag_consensus accepts the per-detector columns and
  ## returns the object with `is_outlier` overwritten under the
  ## requested mode.
  modes <- c("class_aware", "strict", "majority", "speed_trusted", "any")
  cc <- lapply(modes, function(m) {
    n <- tryCatch(
      sum(as.logical(
        move2utils::mt_flag_consensus(x, mode = m)$is_outlier),
        na.rm = TRUE),
      error = function(e) NA_integer_)
    data.frame(consensus_mode = m,
               n_flagged      = n,
               pct_total      = if (is.na(n)) NA_real_
                                else round(100 * n / nrow(ev), 3),
               delta_vs_current = if (is.na(n)) NA_integer_
                                  else n - n_flag,
               stringsAsFactors = FALSE)
  })
  table_consensus_comparison <- do.call(rbind, cc)
  attr(table_consensus_comparison, "current_n_flagged") <- n_flag

  if (print_tables) {
    cat("\n=== Flag audit ===\n")
    cat("Events:", nrow(ev),
        " | flagged (current is_outlier):", n_flag,
        sprintf(" (%.2f%%)\n", 100 * n_flag / max(nrow(ev), 1)))

    cat("\n-- 1. error_class breakdown (flagged only) --\n")
    print(table_error_class, row.names = FALSE)

    cat("\n-- 2. per-detector fires among flagged --\n")
    print(table_detector_fires, row.names = FALSE)

    cat("\n-- 3. co-fire histogram (n detectors fired x is_outlier) --\n")
    print(table_co_fire, row.names = FALSE)

    cat("\n-- 4. near-miss zone (fired but consensus did not trip) --\n")
    print(table_near_miss, row.names = FALSE)

    cat("\n-- 5. consensus mode comparison --\n")
    cat("    Re-applying each mode to existing flag columns;\n")
    cat("    no per-fix detectors are re-run.\n\n")
    print(table_consensus_comparison, row.names = FALSE)
  }

  ## ---- Map: single AEQD projection ----------------------------------
  aeqd <- tryCatch(move2::mt_aeqd_crs(x), error = function(e) NULL)
  x_proj <- if (!is.null(aeqd)) sf::st_transform(x, aeqd) else x
  cc_coords <- sf::st_coordinates(x_proj)

  map_df <- data.frame(
    x = cc_coords[, 1],
    y = cc_coords[, 2],
    tid = as.character(move2::mt_track_id(x_proj)),
    is_outlier = is_out
  )
  map_df$error_class <- if ("error_class" %in% names(ev))
                          as.character(ev$error_class) else NA_character_
  for (col in fire_cols) {
    flagged_v <- !is.na(ev[[col]]) & ev[[col]]
    map_df[[paste0(col, "_status")]] <- ifelse(
      flagged_v & is_out,  "confirmed",
      ifelse(flagged_v & !is_out, "solo", "no_fire"))
  }

  ## Helper to build per-detector panel: kept fixes pale grey, fires
  ## of THIS detector coloured (confirmed = solid, solo = hollow).
  build_detector_panel <- function(col) {
    status_col <- paste0(col, "_status")
    keep_bg <- map_df[map_df[[status_col]] == "no_fire", ]
    fires   <- map_df[map_df[[status_col]] != "no_fire", ]
    label <- sub("^flagged_by_", "", col)
    n_conf <- sum(map_df[[status_col]] == "confirmed")
    n_solo <- sum(map_df[[status_col]] == "solo")
    ggplot2::ggplot() +
      ggplot2::geom_point(data = keep_bg,
                          ggplot2::aes(x = x, y = y),
                          colour = "grey88", size = 0.25, alpha = 0.5) +
      ggplot2::geom_point(data = fires,
                          ggplot2::aes(x = x, y = y, shape = .data[[status_col]]),
                          colour = "#D55E00", size = 1.2, stroke = 0.5,
                          alpha = 0.85) +
      ggplot2::scale_shape_manual(values = c(confirmed = 16, solo = 1),
                                   guide = "none") +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     axis.title       = ggplot2::element_blank(),
                     axis.text        = ggplot2::element_text(size = 5),
                     plot.title       = ggplot2::element_text(size = 8.5)) +
      ggplot2::labs(title = sprintf("%s  (confirmed: %d, solo: %d)",
                                     label, n_conf, n_solo))
  }

  detector_panels <- lapply(fire_cols, build_detector_panel)

  ## Consensus / error_class summary panel (bottom row).
  flagged_pts <- map_df[map_df$is_outlier, ]
  kept_pts    <- map_df[!map_df$is_outlier, ]
  p_consensus <- ggplot2::ggplot() +
    ggplot2::geom_point(data = kept_pts,
                        ggplot2::aes(x = x, y = y),
                        colour = "grey85", size = 0.25, alpha = 0.5) +
    ggplot2::geom_point(data = flagged_pts,
                        ggplot2::aes(x = x, y = y, colour = error_class),
                        size = 1.4, alpha = 0.9) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   axis.title       = ggplot2::element_blank(),
                   axis.text        = ggplot2::element_text(size = 6),
                   plot.title       = ggplot2::element_text(size = 10),
                   legend.position  = "right") +
    ggplot2::labs(title = sprintf(
      "Consensus result (is_outlier): %d flagged (%.2f%%)",
      n_flag, 100 * n_flag / max(nrow(ev), 1)),
      colour = "error_class")

  map <- (patchwork::wrap_plots(detector_panels, nrow = 1)) /
         p_consensus +
    patchwork::plot_layout(heights = c(1, 2)) +
    patchwork::plot_annotation(
      caption = paste(
        "Top row: per-detector fires; solid = confirmed by consensus, hollow = solo fire (near-miss zone).",
        "Bottom: consensus decision coloured by error_class.",
        "Single AEQD projection centred on the move2 object.",
        sep = "  "
      ),
      theme = ggplot2::theme(
        plot.caption = ggplot2::element_text(size = 7, colour = "grey45",
                                              hjust = 0)))

  invisible(list(
    error_class           = table_error_class,
    detector_fires        = table_detector_fires,
    co_fire               = table_co_fire,
    near_miss             = table_near_miss,
    consensus_comparison  = table_consensus_comparison,
    map                   = map
  ))
}

## ggplot's tidy-eval aesthetics reference these names symbolically;
## declare them as globals so R CMD check stays quiet.
utils::globalVariables(c("x", "y", "error_class", ".data"))

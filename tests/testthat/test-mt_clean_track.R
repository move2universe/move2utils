## Tests for mt_clean_track() -- the unified orchestrator.

suppressPackageStartupMessages({
  library(move2); library(sf)
})

read_synthetic <- function() {
  path <- system.file("extdata", "synthetic_tracks.csv.gz",
                       package = "move2utils")
  if (nchar(path) == 0)
    path <- "inst/extdata/synthetic_tracks.csv.gz"
  d <- read.csv(gzfile(path), stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  m <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier",
    crs = 4326)
  m <- m[!sf::st_is_empty(m), ]
  m <- dplyr::arrange(m, move2::mt_track_id(m), move2::mt_time(m))
  m
}

read_ground_truth <- function() {
  path <- system.file("extdata", "synthetic_ground_truth.rds",
                       package = "move2utils")
  if (nchar(path) == 0)
    path <- "inst/extdata/synthetic_ground_truth.rds"
  readRDS(path)
}


test_that("rejects non-move2 input", {
  expect_error(mt_clean_track(data.frame(x = 1)),
               "must be a move2 object")
})

test_that("rejects bad v_max", {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  expect_error(mt_clean_track(m_C, v_max = -1), "positive")
  expect_error(mt_clean_track(m_C, v_max = 0),  "positive")
})

test_that("rejects bad iterations", {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  expect_error(mt_clean_track(m_C, iterations = 0),  "positive")
  expect_error(mt_clean_track(m_C, iterations = -5), "positive")
  expect_error(mt_clean_track(m_C, iterations = "bogus"), "until_clean")
})

test_that("iterations = \"until_clean\" is equivalent to Inf", {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  r_str <- suppressMessages(
    mt_clean_track(m_C, iterations = "until_clean", plot = FALSE, remove = FALSE))
  r_inf <- suppressMessages(
    mt_clean_track(m_C, iterations = Inf, plot = FALSE, remove = FALSE))
  expect_equal(which(r_str$is_outlier), which(r_inf$is_outlier))
})

test_that("preserves original CRS after auto-projection", {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  expect_true(sf::st_is_longlat(m_C))
  res <- suppressMessages(mt_clean_track(m_C, plot = FALSE, remove = FALSE))
  expect_identical(sf::st_crs(res), sf::st_crs(m_C))
  expect_equal(nrow(res), nrow(m_C))
})

test_that("outputs the expected diagnostic columns", {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  res <- suppressMessages(mt_clean_track(m_C, plot = FALSE, remove = FALSE))
  for (col in c("is_outlier", "flagged_by_bridge", "flagged_by_prob",
                "flagged_by_speed", "flag_iteration", "block_id",
                "error_class")) {
    expect_true(col %in% names(res), info = col)
  }

  ## allometric route: mass + mode is mutually exclusive with v_max
  expect_error(suppressMessages(mt_clean_track(m_C, v_max = 30,
                                                mass = 5, mode = "flying",
                                                plot = FALSE)),
               "either")
  expect_error(suppressMessages(mt_clean_track(m_C, mass = 5,
                                                plot = FALSE)),
               "both")
  ## allometric route: derives v_max internally and runs to completion
  res_allo <- suppressMessages(mt_clean_track(m_C, mass = 1, mode = "flying",
                                               plot = FALSE, remove = FALSE))
  expect_true("error_class" %in% names(res_allo))
  expect_type(res$is_outlier,        "logical")
  expect_type(res$flagged_by_bridge, "logical")
  expect_type(res$flag_iteration,    "integer")
  expect_type(res$error_class,       "character")
  ## error_class is NA exactly where is_outlier is FALSE
  expect_identical(is.na(res$error_class), !res$is_outlier)
  ## non-NA labels are from the documented vocabulary
  expect_true(all(stats::na.omit(res$error_class) %in%
                    c("physiological", "block", "consensus",
                      "geometric_spike", "state_anomaly",
                      "kinematic_confluence",
                      "state_transition_buffered")))
})

test_that("flags nothing on clean CPF_B", {
  m <- read_synthetic()
  m_B <- m[move2::mt_track_id(m) == "CPF_B", ]
  res <- suppressMessages(mt_clean_track(m_B, plot = FALSE, remove = FALSE))
  ## conjunction rule on clean data should produce zero flags
  expect_equal(sum(res$is_outlier), 0)
})

test_that("recovers ground truth on CPF_C with zero FPs", {
  m  <- read_synthetic()
  gt <- read_ground_truth()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  res <- suppressMessages(mt_clean_track(m_C, plot = FALSE, remove = FALSE))
  truth <- gt$CPF_C$index
  tp <- sum(res$is_outlier[truth])
  fp <- sum(res$is_outlier) - tp
  ## Conjunction is strict by design; we expect most truth fixes to
  ## still be caught because they are spatially extreme (bridge AND
  ## prob both agree). Zero FP is the important property.
  expect_equal(fp, 0)
  expect_gt(tp, 0)
})

test_that("convergence flag is set appropriately", {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  res <- suppressMessages(mt_clean_track(m_C, plot = FALSE, remove = FALSE))
  expect_true(attr(res, "convergence") %in%
              c("no_new_flags", "max_iterations",
                "active_set_too_small", "flag_fraction_exceeded"))
})

test_that("remove = TRUE drops flagged rows", {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  res_keep <- suppressMessages(mt_clean_track(m_C, plot = FALSE, remove = FALSE))
  res_drop <- suppressMessages(mt_clean_track(m_C, plot = FALSE, remove = TRUE))
  expect_equal(nrow(res_drop),
               nrow(res_keep) - sum(res_keep$is_outlier))
  expect_false(any(res_drop$is_outlier))
})

test_that("multi-individual input is processed per-track", {
  m   <- read_synthetic()
  res <- suppressMessages(mt_clean_track(m, plot = FALSE, remove = FALSE))
  expect_equal(nrow(res), nrow(m))
  expect_setequal(unique(as.character(move2::mt_track_id(res))),
                  c("CPF_A", "CPF_B", "CPF_C", "CPF_D", "CPF_E", "CPF_F"))
})

test_that("user-supplied v_max is honoured", {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  res <- suppressMessages(mt_clean_track(m_C, v_max = 50, plot = FALSE, remove = FALSE))
  expect_equal(attr(res, "v_max_used"), 50)
})

test_that("consensus argument is validated", {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  expect_error(mt_clean_track(m_C, consensus = "bogus", plot = FALSE, remove = FALSE),
               "should be one of")
})

test_that("consensus = majority flags >= strict on CPF_C", {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  res_strict <- suppressMessages(
    mt_clean_track(m_C, consensus = "strict", plot = FALSE, remove = FALSE))
  res_maj <- suppressMessages(
    mt_clean_track(m_C, consensus = "majority", plot = FALSE, remove = FALSE))
  ## majority must catch >= strict (looser rule)
  expect_gte(sum(res_maj$is_outlier), sum(res_strict$is_outlier))
})

test_that("consensus = speed_trusted enables speed-only flags", {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  res <- suppressMessages(
    mt_clean_track(m_C, consensus = "speed_trusted", plot = FALSE, remove = FALSE))
  ## must not error; some fixes should be flagged on a bimodal track
  expect_true(sum(res$is_outlier) > 0)
})

test_that("consensus = any gives the maximum flag count", {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  res_any <- suppressWarnings(suppressMessages(
    mt_clean_track(m_C, consensus = "any", plot = FALSE, remove = FALSE)))
  res_strict <- suppressMessages(
    mt_clean_track(m_C, consensus = "strict", plot = FALSE, remove = FALSE))
  expect_gte(sum(res_any$is_outlier), sum(res_strict$is_outlier))
})

test_that("too-short track degrades gracefully", {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ][1:5, ]
  expect_warning(
    res <- mt_clean_track(m_C, plot = FALSE, remove = FALSE),
    "Too few")
  expect_equal(sum(res$is_outlier), 0)
})

test_that(".entropy_or_dip_gap_threshold_lower prefers entropy when valley exists", {
  set.seed(1)
  x <- c(rnorm(500, mean = 0), rnorm(20, mean = -5))  # clear bimodal
  res <- move2utils:::.entropy_or_dip_gap_threshold_lower(x)
  expect_equal(res$method, "valley")
  expect_true(is.na(res$dip_p))  # entropy was sufficient, no dip test needed
})

test_that(".entropy_or_dip_gap_threshold_lower rejects gap on unimodal", {
  set.seed(1)
  x <- rnorm(500)  # unimodal
  res <- move2utils:::.entropy_or_dip_gap_threshold_lower(x)
  ## Either entropy returned nothing AND gap was rejected by the dip,
  ## or entropy found nothing and gap found nothing -- either way no flags.
  expect_equal(sum(res$is_outlier), 0)
})


# ---- block-expansion gate -------------------------------------------------

test_that(".gate_block_expansion declines when largest does not dominate", {
  ## A trajectory carved into many medium-sized pieces by a too-low
  ## cap.  No single component holds the trajectory; expansion would
  ## flag legitimate movement as blocks.
  sizes <- c(8000L, 6000L, 5000L, 4500L, 3500L, 3000L)
  part  <- list(comp_kept = rep(seq_along(sizes), times = sizes),
                sizes     = sizes,
                kept      = seq_len(sum(sizes)),
                n         = sum(sizes))
  gate <- move2utils:::.gate_block_expansion(part, max_block_fraction = 0.2)
  expect_false(gate$allow)
  expect_match(gate$reason, "severs|trajectory", ignore.case = TRUE)
  expect_lt(gate$largest_frac, 0.8)
})

test_that(".gate_block_expansion declines when sizes form a continuum", {
  ## largest holds 75% (just below the 80% dominance bar) and the rest
  ## form a continuous tail without a clear gap above the bulk -- a
  ## cap that has fragmented sub-trajectories rather than isolated
  ## error clusters.  Either dominance fails or the upper-tail gap
  ## test fails; in either case decline.
  sizes <- as.integer(round(c(75000, 8000, 5000, 3000, 2000, 1500, 1000, 700)))
  part  <- list(comp_kept = rep(seq_along(sizes), times = sizes),
                sizes     = sizes,
                kept      = seq_len(sum(sizes)),
                n         = sum(sizes))
  gate <- move2utils:::.gate_block_expansion(part, max_block_fraction = 0.2)
  expect_false(gate$allow)
})

test_that(".gate_block_expansion allows on trajectory + isolated clusters", {
  ## one giant component (95%) + a handful of tiny isolated clusters
  ## with a clean log-gap above the bulk: the spoof / teleport
  ## signature the gate is designed to permit.
  sizes <- c(50000L, 200L, 180L, 220L, 150L, 250L, 170L, 210L)
  part  <- list(comp_kept = rep(seq_along(sizes), times = sizes),
                sizes     = sizes,
                kept      = seq_len(sum(sizes)),
                n         = sum(sizes))
  gate <- move2utils:::.gate_block_expansion(part, max_block_fraction = 0.2)
  expect_true(gate$allow)
  expect_equal(gate$size_break, 50000L)
  expect_gte(gate$gap_factor, 3)
})

test_that(".gate_block_expansion handles trivial cases", {
  ## single component => nothing to expand
  part1 <- list(comp_kept = rep(1L, 100L), sizes = 100L,
                kept = 1:100, n = 100L)
  g1 <- move2utils:::.gate_block_expansion(part1, max_block_fraction = 0.2)
  expect_false(g1$allow)
  expect_match(g1$reason, "1 component", fixed = TRUE)

  ## two components with one dominant => allow
  part2 <- list(comp_kept = c(rep(1L, 950L), rep(2L, 50L)),
                sizes = c(950L, 50L),
                kept = 1:1000, n = 1000L)
  g2 <- move2utils:::.gate_block_expansion(part2, max_block_fraction = 0.2)
  expect_true(g2$allow)

  ## two components nearly equal in size => decline
  part3 <- list(comp_kept = c(rep(1L, 510L), rep(2L, 490L)),
                sizes = c(510L, 490L),
                kept = 1:1000, n = 1000L)
  g3 <- move2utils:::.gate_block_expansion(part3, max_block_fraction = 0.2)
  expect_false(g3$allow)
})


# ---- silent parameter -----------------------------------------------------

test_that("silent = TRUE suppresses messages from mt_clean_track", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  ## Default (silent = FALSE) emits messages
  expect_message(mt_clean_track(m_A, plot = FALSE, remove = FALSE),
                 "mt_clean_track")
  ## silent = TRUE produces none
  expect_no_message(
    mt_clean_track(m_A, plot = FALSE, remove = FALSE, silent = TRUE))
})

test_that("silent propagates from orchestrator to per-individual dispatch", {
  m <- read_synthetic()
  ## multi-individual call -- silent should suppress the per-id banner too
  expect_no_message(
    mt_clean_track(m, plot = FALSE, remove = FALSE, silent = TRUE))
})


# ---- per-track (mass, mode) ----------------------------------------------

test_that("mass accepts per-track named numeric and dispatches per individual", {
  m <- read_synthetic()
  ## Two tracks: CPF_A and CPF_C.  Pass per-track masses.
  m_AC <- m[move2::mt_track_id(m) %in% c("CPF_A", "CPF_C"), ]
  res <- suppressMessages(mt_clean_track(m_AC,
    mass = c("CPF_A" = 1.5, "CPF_C" = 5.0),
    mode = "flying",
    plot = FALSE, remove = FALSE))
  expect_true("is_outlier" %in% names(res))
  expect_equal(nrow(res), nrow(m_AC))
})

test_that("per-track mass errors when names don't cover all tracks", {
  m <- read_synthetic()
  m_AC <- m[move2::mt_track_id(m) %in% c("CPF_A", "CPF_C"), ]
  expect_error(suppressMessages(mt_clean_track(m_AC,
    mass = c("CPF_A" = 1.5),     # CPF_C missing
    mode = "flying", plot = FALSE)),
    "missing entries for")
})

test_that("mass with a units attribute is auto-converted to kg", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  ## 5 kg expressed as 5000 g; result should match a direct kg call
  mass_g <- units::set_units(5000, "g")
  res_g <- suppressMessages(mt_clean_track(m_A, mass = mass_g, mode = "flying",
                                              plot = FALSE, remove = FALSE))
  res_kg <- suppressMessages(mt_clean_track(m_A, mass = 5,        mode = "flying",
                                              plot = FALSE, remove = FALSE))
  expect_identical(res_g$is_outlier, res_kg$is_outlier)
})

test_that("bare scalar mass > 100 warns about likely-grams", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  expect_warning(
    suppressMessages(mt_clean_track(m_A, mass = 5000, mode = "flying",
                                       plot = FALSE, remove = FALSE)),
    "looks like grams|grams|kg")
})


# ---- compact verbose level ------------------------------------------------

test_that("compact = TRUE suppresses per-iteration output but keeps summary", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  msgs <- testthat::capture_messages(
    mt_clean_track(m_A, plot = FALSE, remove = FALSE, compact = TRUE))
  ## Final summary line is still emitted...
  expect_true(any(grepl("=== mt_clean_track:", msgs)))
  ## ...but no per-iteration "Iter N: bridge=" lines.
  expect_false(any(grepl("^Iter [0-9]+: bridge=", msgs)))
})

test_that("compact = TRUE on multi-individual yields one summary per track", {
  m <- read_synthetic()
  msgs <- testthat::capture_messages(
    mt_clean_track(m, plot = FALSE, remove = FALSE, compact = TRUE))
  n_summary <- sum(grepl("=== mt_clean_track:", msgs))
  n_iter    <- sum(grepl("^Iter [0-9]+: bridge=", msgs))
  ## three tracks => three per-track summary lines
  expect_gte(n_summary, 3L)
  expect_equal(n_iter, 0L)
})

## ---- pre_peel_aux = "primitives" (asymmetric pre-peel) ------------
## Empirical anchor: audit 2026-05-11 (CASCADE_AUDIT_2026-05-11.md
## Section 6.1) showed asymmetric pre-peel rescues the mass-mode CPF_A
## regression and cuts CPF_D's block-boundary FPs.  These tests lock in
## the per-track behaviour so future cascade refactors can't silently
## undo the wins.

test_that("pre_peel_aux = 'none' default reproduces baseline mass-mode behaviour", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  out_default <- suppressMessages(suppressWarnings(
    mt_clean_track(m_A, mass = 1.0, mode = "flying",
                   plot = FALSE, remove = FALSE, silent = TRUE)))
  out_none    <- suppressMessages(suppressWarnings(
    mt_clean_track(m_A, mass = 1.0, mode = "flying", pre_peel_aux = "none",
                   plot = FALSE, remove = FALSE, silent = TRUE)))
  ## Default and explicit "none" must produce identical flag sets.
  expect_identical(out_default$is_outlier, out_none$is_outlier)
})

test_that("pre_peel_aux = 'primitives' improves CPF_A F1 over symmetric default", {
  m <- read_synthetic()
  gt <- read_ground_truth()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  truth <- logical(nrow(m_A))
  truth[as.integer(gt$CPF_A$index)] <- TRUE
  f1 <- function(o) {
    tp <- sum(o & truth); fp <- sum(o & !truth); fn <- sum(!o & truth)
    p <- tp / (tp + fp); r <- tp / (tp + fn)
    if (p + r > 0) 2 * p * r / (p + r) else 0
  }
  out_sym <- suppressMessages(suppressWarnings(
    mt_clean_track(m_A, mass = 1.0, mode = "flying", pre_peel_aux = "none",
                   plot = FALSE, remove = FALSE, silent = TRUE)))
  out_asy <- suppressMessages(suppressWarnings(
    mt_clean_track(m_A, mass = 1.0, mode = "flying", pre_peel_aux = "primitives",
                   plot = FALSE, remove = FALSE, silent = TRUE)))
  ## Audit anchor: CPF_A symmetric F1 = 0.885; asymmetric F1 = 0.958.
  ## Test margin (>= 0.04) absorbs minor cascade tweaks while still
  ## guarding the round-3 #1 win.
  expect_gt(f1(out_asy$is_outlier), f1(out_sym$is_outlier) + 0.04)
})

test_that("pre_peel_aux = 'primitives' improves CPF_D F1 over symmetric default", {
  m <- read_synthetic()
  gt <- read_ground_truth()
  m_D <- m[move2::mt_track_id(m) == "CPF_D", ]
  truth <- logical(nrow(m_D))
  truth[as.integer(gt$CPF_D$index)] <- TRUE
  f1 <- function(o) {
    tp <- sum(o & truth); fp <- sum(o & !truth); fn <- sum(!o & truth)
    p <- tp / (tp + fp); r <- tp / (tp + fn)
    if (p + r > 0) 2 * p * r / (p + r) else 0
  }
  out_sym <- suppressMessages(suppressWarnings(
    mt_clean_track(m_D, mass = 1.0, mode = "flying", pre_peel_aux = "none",
                   plot = FALSE, remove = FALSE, silent = TRUE)))
  out_asy <- suppressMessages(suppressWarnings(
    mt_clean_track(m_D, mass = 1.0, mode = "flying", pre_peel_aux = "primitives",
                   plot = FALSE, remove = FALSE, silent = TRUE)))
  ## Audit anchor: CPF_D symmetric F1 = 0.674; asymmetric F1 ~ 0.78.
  expect_gt(f1(out_asy$is_outlier), f1(out_sym$is_outlier) + 0.06)
})

## ---- block-expansion monotonicity in user information --------------
## Regression for the seam-dilution / `!used_peel` bug (2026-06-23):
## supplying a physiological cap must NOT reduce recovery of a coherent
## boundary block.  Before the fix the auto path recovered 150/150 (block
## expansion) while `mass`/`mode` recovered 2/150 -- the peel removed the
## seam, widened the gap, and the diluted across-seam speed fell below the
## cap, so block expansion found nothing.  The fix computes the partition
## on original (un-diluted) timing and drops the `!used_peel` short-circuit.
## Fixture = sustained boundary-anchored spoof block
## (inst/extdata/make_boundary_spoof_demo.R); this is the ONLY construction
## that exercises graph block-expansion (CPF_D's mid-track block is
## recovered by the per-fix layer, with the gate correctly declining).

test_that("supplying a physiological cap does not reduce block recovery", {
  src <- system.file("extdata", "make_boundary_spoof_demo.R",
                     package = "move2utils")
  if (nchar(src) == 0) src <- "inst/extdata/make_boundary_spoof_demo.R"
  source(src, local = TRUE)
  d <- make_boundary_spoof_demo()
  nb <- length(d$truth)

  o_auto <- suppressMessages(suppressWarnings(
    mt_clean_track(d$track, plot = FALSE, remove = FALSE, silent = TRUE)))
  o_cap  <- suppressMessages(suppressWarnings(
    mt_clean_track(d$track, mass = 1, mode = "flying",
                   plot = FALSE, remove = FALSE, silent = TRUE)))

  rec_auto <- sum(o_auto$is_outlier[d$truth])
  rec_cap  <- sum(o_cap$is_outlier[d$truth])

  ## auto path recovers essentially the whole block via block-expansion
  expect_gt(rec_auto, 0.9 * nb)
  ## monotonicity: the cap path must do at least as well (was 2/nb)
  expect_gte(rec_cap, rec_auto)
  ## and it is genuinely block-expansion doing the work on the cap path
  expect_gt(sum(!is.na(o_cap$block_id)), 0.9 * nb)
})

test_that("block expansion declines on a clean track under a physiological cap", {
  ## The spike-fragmentation caveat: a clean trajectory under a cap must
  ## not manufacture blocks (no minority component dominates).
  m <- read_synthetic()
  m_B <- m[move2::mt_track_id(m) == "CPF_B", ]
  o <- suppressMessages(suppressWarnings(
    mt_clean_track(m_B, mass = 1, mode = "flying",
                   plot = FALSE, remove = FALSE, silent = TRUE)))
  expect_identical(sum(!is.na(o$block_id)), 0L)
})

## ---- persistence_filter = "class_aware" (R2 opt-in) ---------------
## Empirical anchor (CASCADE_AUDIT_2026-05-11.md Section 6.2 + the
## 2026-05-09 class-conditional finding, post-CRS-fix numbers):
## the filter targets only `state_anomaly` + `consensus` classes and
## demotes flags whose persistence_count < 3.  Net CPF effect is
## near-neutral (mean DF1 = -0.003); shipping as opt-in only.

test_that("persistence_filter = 'none' default is identity", {
  m <- read_synthetic()
  m_D <- m[move2::mt_track_id(m) == "CPF_D", ]
  out_default <- suppressMessages(suppressWarnings(
    mt_clean_track(m_D, plot = FALSE, remove = FALSE, silent = TRUE)))
  out_none    <- suppressMessages(suppressWarnings(
    mt_clean_track(m_D, persistence_filter = "none",
                   plot = FALSE, remove = FALSE, silent = TRUE)))
  expect_identical(out_default$is_outlier, out_none$is_outlier)
})

test_that("persistence_filter = 'class_aware' demotes only state_anomaly + consensus", {
  m <- read_synthetic()
  m_D <- m[move2::mt_track_id(m) == "CPF_D", ]
  out_n <- suppressMessages(suppressWarnings(
    mt_clean_track(m_D, plot = FALSE, remove = FALSE, silent = TRUE)))
  out_p <- suppressMessages(suppressWarnings(
    mt_clean_track(m_D, persistence_filter = "class_aware",
                   plot = FALSE, remove = FALSE, silent = TRUE)))
  ## Filter only ever DROPS flags; never adds.
  expect_true(all(out_p$is_outlier <= out_n$is_outlier))
  ## Any dropped flag was in state_anomaly OR consensus class
  ## (relative to the unfiltered output's classification).
  dropped <- which(out_n$is_outlier & !out_p$is_outlier)
  if (length(dropped))
    expect_true(all(out_n$error_class[dropped] %in%
                     c("state_anomaly", "consensus")))
})

## ---- max_iterations default = 100L (R8 resolution) -----------------
## Per CASCADE_AUDIT_2026-05-11 Section 3.4 + 8: the iteration-count
## distribution across all audit tracks tops out at 22 iter for
## legitimate convergence; K02 (block contamination) is the
## longest at 61 iter.  Default 100 covers K02 with margin and caps
## multi-state non-convergence (WH17-class) at half the previous
## wallclock.

test_that("default max_iterations is 100", {
  expect_equal(formals(mt_clean_track)$max_iterations, 100L)
})

test_that("pre_peel_aux is no-op without v_max / mass-mode", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  ## With no physiological cap the pre-peel block does not run, so
  ## "primitives" cannot change the outcome -- result must equal the
  ## auto-cap baseline.
  out_auto  <- suppressMessages(suppressWarnings(
    mt_clean_track(m_A, plot = FALSE, remove = FALSE, silent = TRUE)))
  out_aux   <- suppressMessages(suppressWarnings(
    mt_clean_track(m_A, pre_peel_aux = "primitives",
                   plot = FALSE, remove = FALSE, silent = TRUE)))
  expect_identical(out_auto$is_outlier, out_aux$is_outlier)
})


## ---- pool_by orchestrator tests --------------------------------------

## Helper: 2-track move2 sharing one animal_id, for the orchestrator
## pool_by suite.  Cannot reuse read_synthetic() since CPF tracks are
## different animals.
make_pair_orch <- function(n1 = 200, n2 = 100, indv = "I1") {
  set.seed(11)
  build <- function(id, n, lon0, lat0) {
    dx <- rnorm(n, 0, 0.002); dy <- rnorm(n, 0, 0.002)
    data.frame(id = id,
               timestamp = as.POSIXct("2026-01-01", tz = "UTC") +
                            seq_len(n) * 3600,
               lon = lon0 + cumsum(dx), lat = lat0 + cumsum(dy))
  }
  d <- rbind(build("t1", n1, 10, 50), build("t2", n2, 12, 51))
  m <- move2::mt_as_move2(d, coords = c("lon", "lat"),
                            time_column = "timestamp",
                            track_id_column = "id", crs = 4326)
  td <- move2::mt_track_data(m); td$indv <- indv
  move2::mt_set_track_data(m, td)
}


test_that("orchestrator pool_by validates inputs", {
  m <- make_pair_orch()
  ## non-character / wrong shape: helpful message names the 2-level cap.
  expect_error(mt_clean_track(m, pool_by = 1, plot = FALSE, remove = FALSE),
               "length 1 \\(single column")
  expect_error(mt_clean_track(m, pool_by = c("a","b","c"),
                                plot = FALSE, remove = FALSE),
               "Deeper hierarchies")
  ## non-existent columns: error names them.
  expect_error(suppressMessages(
                  mt_clean_track(m, pool_by = "no_such_col",
                                   plot = FALSE, remove = FALSE, silent = TRUE)),
               "not in")
  expect_error(suppressMessages(
                  mt_clean_track(m, pool_by = c("a","b"),
                                   plot = FALSE, remove = FALSE, silent = TRUE)),
               "not in")
  ## identical outer == inner under length-2 form.
  expect_error(mt_clean_track(m, pool_by = c("indv","indv"),
                                plot = FALSE, remove = FALSE, silent = TRUE),
               "two \\*distinct\\* columns")
})


test_that("orchestrator pool_by = NULL is byte-identical to no pool_by", {
  m <- make_pair_orch()
  o1 <- suppressMessages(suppressWarnings(
          mt_clean_track(m, plot = FALSE, remove = FALSE, silent = TRUE)))
  o2 <- suppressMessages(suppressWarnings(
          mt_clean_track(m, pool_by = NULL, plot = FALSE, remove = FALSE,
                          silent = TRUE)))
  expect_identical(o1$is_outlier, o2$is_outlier)
  expect_identical(o1$flagged_by_bridge, o2$flagged_by_bridge)
  expect_identical(o1$flagged_by_detour, o2$flagged_by_detour)
})


test_that("orchestrator pool_by union is strictly additive", {
  ## per-track flags ⊆ pool flags; pool path is purely additive.
  m <- make_pair_orch()
  o_pt <- suppressMessages(suppressWarnings(
            mt_clean_track(m, plot = FALSE, remove = FALSE, silent = TRUE)))
  o_pg <- suppressMessages(suppressWarnings(
            mt_clean_track(m, pool_by = "indv",
                            plot = FALSE, remove = FALSE, silent = TRUE)))
  expect_true(all(which(o_pt$is_outlier) %in% which(o_pg$is_outlier)))
})


test_that("orchestrator pool-added flags carry valid error_class", {
  ## Pool-added fixes get error_class = "pool" when the cascade
  ## didn't classify them; cascade-classified fixes keep their class.
  m <- make_pair_orch(n1 = 300, n2 = 80)
  o <- suppressMessages(suppressWarnings(
          mt_clean_track(m, pool_by = "indv",
                          plot = FALSE, remove = FALSE, silent = TRUE)))
  expect_true("error_class" %in% names(o))
  ## All flagged fixes must have a known error_class (including "pool").
  known <- c("pool", "consensus", "geometric_spike",
              "state_anomaly", "kinematic_confluence",
              "block", "physiological",
              "state_transition_buffered")
  flagged_classes <- o$error_class[o$is_outlier]
  expect_true(all(is.na(flagged_classes) | flagged_classes %in% known))
})


test_that("orchestrator pool_by + state errors on disjoint state vocab", {
  set.seed(33); n <- 80
  build <- function(id, n, lon0, lat0, state_label) {
    data.frame(id = id,
               timestamp = as.POSIXct("2026-01-01", tz = "UTC") +
                            seq_len(n) * 3600,
               lon = lon0 + rnorm(n, 0, 0.002),
               lat = lat0 + rnorm(n, 0, 0.002),
               state = state_label)
  }
  d <- rbind(build("t1", n, 10, 50, "flight"),
              build("t2", n, 12, 51, "sleep"))
  m <- move2::mt_as_move2(d, coords = c("lon","lat"),
                            time_column = "timestamp",
                            track_id_column = "id", crs = 4326)
  td <- move2::mt_track_data(m); td$indv <- "I1"
  m <- move2::mt_set_track_data(m, td)
  expect_error(mt_clean_track(m, state = "state", pool_by = "indv",
                                plot = FALSE, remove = FALSE,
                                silent = TRUE),
               "vocabulary mismatch")
})


test_that("orchestrator pool_by + state OK when vocabularies intersect", {
  set.seed(33); n <- 80
  build <- function(id, n, lon0, lat0, states) {
    data.frame(id = id,
               timestamp = as.POSIXct("2026-01-01", tz = "UTC") +
                            seq_len(n) * 3600,
               lon = lon0 + rnorm(n, 0, 0.002),
               lat = lat0 + rnorm(n, 0, 0.002),
               state = states)
  }
  d <- rbind(build("t1", n, 10, 50, rep("flight", n)),
              build("t2", n, 12, 51, c(rep("flight", n / 2),
                                        rep("rest",   n / 2))))
  m <- move2::mt_as_move2(d, coords = c("lon","lat"),
                            time_column = "timestamp",
                            track_id_column = "id", crs = 4326)
  td <- move2::mt_track_data(m); td$indv <- "I1"
  m <- move2::mt_set_track_data(m, td)
  expect_no_error(suppressMessages(suppressWarnings(
    mt_clean_track(m, state = "state", pool_by = "indv",
                    plot = FALSE, remove = FALSE, silent = TRUE))))
})


## ---- Item G / R7 primitive-knob overrides ----------------------------

test_that("primitive-knob overrides: all NULL is byte-identical to defaults", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  o_def <- suppressMessages(suppressWarnings(
            mt_clean_track(m_A, plot = FALSE, remove = FALSE, silent = TRUE)))
  o_null <- suppressMessages(suppressWarnings(
            mt_clean_track(m_A, plot = FALSE, remove = FALSE, silent = TRUE,
                            bridge_method         = NULL,
                            bridge_threshold_type = NULL,
                            bridge_iterations     = NULL,
                            prob_threshold_type   = NULL,
                            detour_threshold_type = NULL)))
  expect_identical(o_def$is_outlier, o_null$is_outlier)
  expect_identical(o_def$flagged_by_bridge, o_null$flagged_by_bridge)
  expect_identical(o_def$flagged_by_detour, o_null$flagged_by_detour)
  expect_identical(o_def$flagged_by_prob,   o_null$flagged_by_prob)
})


test_that("primitive-knob overrides: bridge_method produces distinct flag sets", {
  ## Each method has a different scope of validity per the paper;
  ## flag sets are not byte-identical across methods.  This test
  ## confirms the override is wired (not silently ignored).
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  o_combined    <- suppressMessages(suppressWarnings(
    mt_clean_track(m_A, bridge_method = "combined",    plot = FALSE,
                    remove = FALSE, silent = TRUE)))
  o_isotropic   <- suppressMessages(suppressWarnings(
    mt_clean_track(m_A, bridge_method = "isotropic",   plot = FALSE,
                    remove = FALSE, silent = TRUE)))
  o_directional <- suppressMessages(suppressWarnings(
    mt_clean_track(m_A, bridge_method = "directional", plot = FALSE,
                    remove = FALSE, silent = TRUE)))
  ## At least one of the alternative methods must differ from combined.
  diff_iso <- !identical(o_combined$is_outlier, o_isotropic$is_outlier)
  diff_dir <- !identical(o_combined$is_outlier, o_directional$is_outlier)
  expect_true(diff_iso || diff_dir,
               info = "Override is silently ignored if all methods match.")
})


test_that("primitive-knob overrides: validation rejects bogus values", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  expect_error(suppressMessages(suppressWarnings(
    mt_clean_track(m_A, bridge_method = "bogus",
                    plot = FALSE, silent = TRUE))),
               "should be one of")
  expect_error(suppressMessages(suppressWarnings(
    mt_clean_track(m_A, bridge_threshold_type = "bogus",
                    plot = FALSE, silent = TRUE))),
               "should be one of")
  expect_error(suppressMessages(suppressWarnings(
    mt_clean_track(m_A, bridge_iterations = -1,
                    plot = FALSE, silent = TRUE))),
               "positive integer")
  expect_error(suppressMessages(suppressWarnings(
    mt_clean_track(m_A, prob_threshold_type = "bogus",
                    plot = FALSE, silent = TRUE))),
               "should be one of")
  expect_error(suppressMessages(suppressWarnings(
    mt_clean_track(m_A, detour_threshold_type = "bogus",
                    plot = FALSE, silent = TRUE))),
               "should be one of")
})


test_that("primitive-knob overrides: bridge_method='directional' on CPF_A wins F1", {
  ## Empirical confirmation that override has meaningful effect (not
  ## just non-trivial diff).  On CPF_A the directional bridge method
  ## eliminates the cascade's 2 FPs without sacrificing TPs.
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  gt  <- read_ground_truth()
  truth <- as.integer(gt$CPF_A$index)
  truth <- truth[truth >= 1 & truth <= nrow(m_A)]
  is_truth <- logical(nrow(m_A)); is_truth[truth] <- TRUE

  o <- suppressMessages(suppressWarnings(
    mt_clean_track(m_A, bridge_method = "directional",
                    plot = FALSE, remove = FALSE, silent = TRUE)))
  tp <- sum(o$is_outlier & is_truth)
  fp <- sum(o$is_outlier & !is_truth)
  fn <- sum(!o$is_outlier & is_truth)
  prec <- tp / max(tp + fp, 1); rec <- tp / max(tp + fn, 1)
  f1 <- 2 * prec * rec / max(prec + rec, 1e-9)
  expect_equal(tp, 23L)
  expect_lte(fp, 1L)  # directional should produce <= 1 FP (sweep-validated 0)
  expect_gte(f1, 0.97)
})

test_that("combined_evidence is exposed under evidence modes only", {
  m   <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  ec <- suppressMessages(mt_clean_track(m_A, remove = FALSE, plot = FALSE,
                                        silent = TRUE))  # default evidence_corroborated
  ca <- suppressMessages(mt_clean_track(m_A, remove = FALSE, plot = FALSE,
                                        silent = TRUE, consensus = "class_aware"))
  ## evidence modes attach the one-lever score; Boolean modes do not
  expect_true("combined_evidence" %in% names(ec))
  expect_type(ec$combined_evidence, "double")
  expect_false("combined_evidence" %in% names(ca))
  ## flagged fixes carry higher evidence than kept ones
  fl <- ec$is_outlier %in% TRUE
  expect_gt(mean(ec$combined_evidence[fl], na.rm = TRUE),
            mean(ec$combined_evidence[!fl], na.rm = TRUE))
})

## Regression + extension tests for the generalised consensus contract.
##
## The behaviour invariant of the 2026-06-05 refactor: with the four
## standard detectors and no competence weights, every legacy mode is
## byte-identical to the pre-refactor four-vector implementation.  The
## reference implementation below encodes that original logic verbatim;
## the tests assert the new matrix core matches it on random inputs.

## --- reference implementation (the original .consensus_logical) -------
ref_consensus <- function(b, p, s, d, mode) {
  geo   <- b | d
  votes <- as.integer(b) + as.integer(p) + as.integer(s) + as.integer(d)
  switch(mode,
    class_aware   = (votes >= 3L) | (b & d) | (geo & s) | (geo & p),
    strict        = geo & (p | s),
    majority      = votes >= 2L,
    speed_trusted = s | (geo & p),
    any           = b | p | s | d
  )
}

random_flags <- function(n, seed) {
  set.seed(seed)
  list(b = as.logical(rbinom(n, 1, 0.25)),
       p = as.logical(rbinom(n, 1, 0.25)),
       s = as.logical(rbinom(n, 1, 0.25)),
       d = as.logical(rbinom(n, 1, 0.25)))
}

test_that("legacy modes are byte-identical to the reference on random input", {
  for (seed in 1:20) {
    f <- random_flags(200, seed)
    for (mode in c("class_aware", "strict", "majority",
                   "speed_trusted", "any")) {
      got <- move2utils:::.consensus_logical(f$b, f$p, f$s, f$d, mode = mode)
      exp <- ref_consensus(f$b, f$p, f$s, f$d, mode)
      expect_identical(got, exp,
                       info = paste("mode", mode, "seed", seed))
    }
  }
})

test_that(".consensus_decide matches .consensus_logical via the matrix path", {
  f <- random_flags(100, 99)
  flags <- cbind(bridge = f$b, prob = f$p, speed = f$s, detour = f$d)
  for (mode in c("class_aware", "strict", "majority",
                 "speed_trusted", "any")) {
    via_matrix <- move2utils:::.consensus_decide(flags, mode = mode)
    via_legacy <- move2utils:::.consensus_logical(f$b, f$p, f$s, f$d, mode = mode)
    expect_identical(via_matrix, via_legacy, info = mode)
  }
})

test_that("competence_count with NULL weights equals majority", {
  f <- random_flags(150, 7)
  flags <- cbind(bridge = f$b, prob = f$p, speed = f$s, detour = f$d)
  cc  <- move2utils:::.consensus_decide(flags, mode = "competence_count", k = 2L)
  maj <- move2utils:::.consensus_decide(flags, mode = "majority")
  expect_identical(cc, maj)
})

test_that("competence_count silences a detector below the gate", {
  ## Three fixes; only the bridge+detour pair fires, but detour is
  ## incompetent (weight 0) at fix 2, so only fix 1 survives at k = 2.
  flags <- cbind(bridge = c(TRUE,  TRUE,  FALSE),
                 prob   = c(FALSE, FALSE, FALSE),
                 speed  = c(FALSE, FALSE, FALSE),
                 detour = c(TRUE,  TRUE,  FALSE))
  w     <- cbind(bridge = c(1, 1, 1),
                 prob   = c(1, 1, 1),
                 speed  = c(1, 1, 1),
                 detour = c(1, 0, 1))   # detour silenced at fix 2
  out <- move2utils:::.consensus_decide(flags, weights = w,
                                        mode = "competence_count",
                                        gate = 0.5, k = 2L)
  expect_identical(out, c(TRUE, FALSE, FALSE))
})

test_that("vote-counting modes generalise to a fifth detector", {
  ## A fifth detector adds to the vote count for majority/any.
  flags <- cbind(bridge = c(TRUE,  FALSE),
                 prob   = c(TRUE,  FALSE),
                 speed  = c(FALSE, FALSE),
                 detour = c(FALSE, FALSE),
                 speed_rel = c(FALSE, TRUE))
  ## fix 1: 2 votes -> majority TRUE; fix 2: 1 vote -> majority FALSE, any TRUE
  expect_identical(move2utils:::.consensus_decide(flags, mode = "majority"),
                   c(TRUE, FALSE))
  expect_identical(move2utils:::.consensus_decide(flags, mode = "any"),
                   c(TRUE, TRUE))
})

test_that("legacy 4-arg custom is still honoured by the core", {
  f <- random_flags(50, 3)
  flags <- cbind(bridge = f$b, prob = f$p, speed = f$s, detour = f$d)
  out <- move2utils:::.consensus_decide(
    flags, mode = "custom", custom = function(b, p, s, d) b & s)
  expect_identical(out, f$b & f$s)
})

test_that("general 2-arg custom receives the flag matrix", {
  flags <- cbind(bridge = c(TRUE, FALSE), prob = c(FALSE, TRUE),
                 speed = c(FALSE, FALSE), detour = c(FALSE, FALSE))
  out <- move2utils:::.consensus_decide(
    flags, mode = "custom",
    custom = function(flags, weights) rowSums(flags) >= 1L)
  expect_identical(out, c(TRUE, TRUE))
})

test_that("weighted_evidence sums log-LRs and thresholds the posterior", {
  ## three detectors, three fixes
  ev <- cbind(bridge = c( 2.0, -1.0,  0.0),
              detour = c( 0.5,  0.5,  3.0),
              prob   = c(-0.1, -0.1, -0.1))
  ## row sums: 2.4, -0.6, 2.9  -> >0 at fix 1 and 3
  out <- move2utils:::.consensus_decide(
    flags = matrix(FALSE, 3, 0), mode = "weighted_evidence", evidence = ev)
  expect_identical(out, c(TRUE, FALSE, TRUE))
})

test_that("weighted_evidence treats NA as abstention (zero), not evidence-against", {
  ## abstaining bridge (NA) must not drag the posterior down
  ev_abstain <- cbind(bridge = c(NA, NA), detour = c(1.5, -0.2))
  out <- move2utils:::.consensus_decide(
    flags = matrix(FALSE, 2, 0), mode = "weighted_evidence", evidence = ev_abstain)
  expect_identical(out, c(TRUE, FALSE))   # fix1: 0+1.5>0 ; fix2: 0-0.2<0
  ## contrast: a competent silent bridge (negative log-LR) CAN veto
  ev_veto <- cbind(bridge = c(-2.0, -2.0), detour = c(1.5, -0.2))
  out2 <- move2utils:::.consensus_decide(
    flags = matrix(FALSE, 2, 0), mode = "weighted_evidence", evidence = ev_veto)
  expect_identical(out2, c(FALSE, FALSE)) # -0.5 and -2.2 both < 0
})

test_that("weighted_evidence honours a non-zero threshold", {
  ev <- cbind(a = c(1.0, 3.0), b = c(0.5, 0.5))
  out <- move2utils:::.consensus_decide(
    flags = matrix(FALSE, 2, 0), mode = "weighted_evidence",
    evidence = ev, evidence_threshold = 2)
  expect_identical(out, c(FALSE, TRUE))   # 1.5 < 2 ; 3.5 > 2
})

test_that("weighted_evidence errors when no evidence supplied", {
  expect_error(
    move2utils:::.consensus_decide(flags = matrix(FALSE, 2, 0),
                                   mode = "weighted_evidence"),
    class = "move2utils_mt_flag_consensus_missing_evidence")
})

test_that("class_aware default is unaffected by the new evidence params", {
  ## adding evidence machinery must not perturb the default decision
  f <- random_flags(80, 5)
  flags <- cbind(bridge = f$b, prob = f$p, speed = f$s, detour = f$d)
  base <- move2utils:::.consensus_decide(flags, mode = "class_aware")
  withev <- move2utils:::.consensus_decide(
    flags, mode = "class_aware",
    evidence = matrix(99, nrow(flags), 1), evidence_threshold = -5)
  expect_identical(base, withev)
})

test_that(".combine_evidence calibration bounds every detector to +/- C", {
  ## one tiny-scale detector, one huge-scale: after calibration both are
  ## bounded to +/- C, so neither dominates.
  ev <- cbind(small = c(-0.5, 0.2, 0.4, -0.1),
              huge  = c(-9e6, 1e3, 8e8, -4e3))
  raw  <- move2utils:::.combine_evidence(ev, calibrate = FALSE)
  cal  <- move2utils:::.combine_evidence(ev, calibrate = TRUE, C = 4)
  ## raw sum is utterly dominated by `huge`; calibrated is bounded to +/-8
  expect_true(max(abs(raw)) > 1e6)
  expect_true(all(abs(cal) <= 8 + 1e-9))
})

test_that(".combine_evidence calibrates per group, never across groups", {
  ## two groups with very different scales; per-group calibration must
  ## scale each group by its OWN MAD, not a pooled global MAD.
  ev    <- cbind(a = c(1, 2, 3,  10, 20, 30))   # group B is 10x group A
  group <- c("A","A","A","B","B","B")
  gcal  <- move2utils:::.combine_evidence(ev, calibrate = TRUE, C = 4, group = group)
  ## computing each group alone must match the grouped result on its rows
  aonly <- move2utils:::.combine_evidence(ev[1:3, , drop = FALSE], calibrate = TRUE, C = 4)
  bonly <- move2utils:::.combine_evidence(ev[4:6, , drop = FALSE], calibrate = TRUE, C = 4)
  expect_equal(gcal[1:3], aonly)
  expect_equal(gcal[4:6], bonly)
  ## and a global (ungrouped) calibration would differ from the per-group one
  glob <- move2utils:::.combine_evidence(ev, calibrate = TRUE, C = 4)
  expect_false(isTRUE(all.equal(glob, gcal)))
})

test_that(".loglr_grouped computes the flag boundary per group", {
  np  <- c(1, 2, 5,   1, 2, 9)      # surprisal
  io  <- c(F, F, T,   F, F, T)      # one outlier per group
  grp <- c("A","A","A","B","B","B")
  g  <- move2utils:::.loglr_grouped(np, io, grp)
  ## per group the boundary is the flagged fix's surprisal -> its loglr 0
  expect_equal(g[3], 0)             # group A boundary at 5
  expect_equal(g[6], 0)             # group B boundary at 9 (NOT the global 5)
})

test_that(".combine_evidence treats NA as abstention under calibration", {
  ev <- cbind(a = c(NA, 2, -2), b = c(1.5, NA, NA))
  cal <- move2utils:::.combine_evidence(ev, calibrate = TRUE, C = 4)
  expect_length(cal, 3)
  expect_true(all(is.finite(cal)))     # NA -> 0 contribution, never NA out
})

test_that("weighted_evidence + calibrate attaches combined_evidence and decides on it", {
  suppressPackageStartupMessages({library(move2); library(sf)})
  n <- 10
  d <- data.frame(t = as.POSIXct("2025-01-01", tz = "UTC") + seq_len(n) * 3600,
                  id = "a")
  g <- sf::st_sfc(lapply(seq_len(n), function(i) sf::st_point(c(i, i))), crs = 4326)
  x <- mt_as_move2(sf::st_sf(d, geometry = g), time_column = "t",
                   track_id_column = "id")
  ## two log-LR columns on wildly different scales
  x$loglr_a <- c(-1, -1, -1, 5e6, -1, 3, -1, -1, -1, -1)
  x$loglr_b <- c(-0.2, -0.2, -0.2, 0.1, -0.2, 0.5, 0.3, -0.2, -0.2, -0.2)
  out <- mt_flag_consensus(x, mode = "weighted_evidence",
                           evidence_cols = c(a = "loglr_a", b = "loglr_b"),
                           calibrate = TRUE)
  expect_true("combined_evidence" %in% names(out))
  expect_length(out$combined_evidence, n)
  expect_true(out$is_outlier[4]); expect_true(out$is_outlier[6])
  expect_false(out$is_outlier[1])
  ## calibration must stop the 5e6 from making the combined evidence huge
  expect_true(max(abs(out$combined_evidence)) <= 8 + 1e-9)
})

test_that("custom returning wrong length errors with a classed condition", {
  flags <- cbind(bridge = c(TRUE, FALSE), prob = c(FALSE, TRUE),
                 speed = c(FALSE, FALSE), detour = c(FALSE, FALSE))
  expect_error(
    move2utils:::.consensus_decide(flags, mode = "custom",
                                   custom = function(b, p, s, d) TRUE),
    class = "move2utils_mt_flag_consensus_bad_custom_return")
})

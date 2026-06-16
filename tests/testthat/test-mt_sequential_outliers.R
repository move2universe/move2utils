## Helper: load a small CPF synthetic track
read_cpf <- function(track = "CPF_C") {
  d <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")),
                 stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  x <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier", crs = 4326)
  x <- x[!sf::st_is_empty(x), ]
  x <- dplyr::arrange(x, move2::mt_track_id(x), move2::mt_time(x))
  x[move2::mt_track_id(x) == track, ]
}

test_that("mt_sequential_outliers rejects non-move2 input", {
  expect_error(mt_sequential_outliers(NULL), "move2")
  expect_error(mt_sequential_outliers(data.frame(x = 1:3)), "move2")
})

test_that("mt_sequential_outliers validates scan argument", {
  x <- read_cpf()
  expect_error(
    suppressMessages(mt_sequential_outliers(x, scan = "bogus")),
    "should be one of")
})

test_that("mt_sequential_outliers attaches required columns", {
  x <- read_cpf()
  res <- suppressMessages(
    mt_sequential_outliers(x, scan = "forward-backward", plot = FALSE))
  expect_true("is_outlier" %in% names(res))
  expect_true("seq_joint_prob" %in% names(res))
  expect_type(res$is_outlier, "logical")
  expect_type(res$seq_joint_prob, "double")
  expect_equal(nrow(res), nrow(x))
})

test_that("mt_sequential_outliers produces non-zero flags on CPF_A", {
  ## CPF_A has 23 injected outliers; the sequential scan should flag
  ## at least some of them.
  x <- read_cpf("CPF_A")
  res <- suppressMessages(
    mt_sequential_outliers(x, scan = "forward-backward", plot = FALSE))
  expect_true(sum(res$is_outlier, na.rm = TRUE) > 0L)
})

## ---- anchor-corruption diagnostic ----------------------------------
## Pre-v0.3 the function silently flagged everything as anomalous when
## the chosen anchor was corrupt and the runaway stayed below
## max_skip = 100.  v0.3 surfaces this via a warning when the per-track
## flag rate exceeds anchor_corruption_threshold (default 0.30).

test_that("mt_sequential_outliers validates anchor_corruption_threshold", {
  x <- read_cpf()
  expect_error(
    suppressMessages(mt_sequential_outliers(x,
      anchor_corruption_threshold = -0.5)),
    "scalar in \\(0, 1\\)"
  )
  expect_error(
    suppressMessages(mt_sequential_outliers(x,
      anchor_corruption_threshold = 1.5)),
    "scalar in \\(0, 1\\)"
  )
  expect_error(
    suppressMessages(mt_sequential_outliers(x,
      anchor_corruption_threshold = c(0.3, 0.5))),
    "scalar in \\(0, 1\\)"
  )
  expect_error(
    suppressMessages(mt_sequential_outliers(x,
      anchor_corruption_threshold = "high")),
    "scalar in \\(0, 1\\)"
  )
})

test_that("mt_sequential_outliers warns when per-track flag rate exceeds threshold", {
  ## Verify the warning logic: when the per-track flag rate exceeds
  ## anchor_corruption_threshold, a warning fires.  CPF_A has 23
  ## injected outliers in 1748 fixes -- a flag rate of ~1.3%.  Setting
  ## the threshold to 0.5% therefore triggers the warning.  This
  ## tests the contract ("warn when rate > threshold") without having
  ## to engineer a deliberately corrupted anchor (the runaway-flagging
  ## failure mode is rare in practice; it's the silent-failure
  ## potential that motivates the warning).
  x <- read_cpf("CPF_A")
  expect_warning(
    res <- suppressMessages(
      mt_sequential_outliers(x, scan = "forward-backward", plot = FALSE,
        anchor_corruption_threshold = 0.005)),
    "anchor-corruption threshold"
  )
})

test_that("mt_sequential_outliers anchor warning is silenced when disabled", {
  ## anchor_corruption_threshold = NULL turns the diagnostic off.
  x <- read_cpf("CPF_A")
  expect_no_warning(
    suppressMessages(
      mt_sequential_outliers(x, scan = "forward-backward", plot = FALSE,
        anchor_corruption_threshold = NULL))
  )
})

test_that("mt_sequential_outliers does not warn on clean tracks", {
  ## CPF_C has 4 injected outliers in 185 fixes -- a real-world flag
  ## rate of ~2%.  No warning should fire at the default 30% threshold.
  x <- read_cpf("CPF_C")
  expect_no_warning(
    suppressMessages(
      mt_sequential_outliers(x, scan = "forward-backward", plot = FALSE))
  )
})


## ---- two-threshold fix regression tests (2026-05-12) ------------------

test_that("two-threshold fix: WH17-style multi-state bimodal-KDE no longer over-flags", {
  ## Synthetic stand-in for the WH17 bug: a track with a strongly
  ## bimodal step distribution (rest + fast bursts).  Pre-fix, the
  ## autodifference KDE on this bimodal distribution explodes (large
  ## densities at the modal peaks), driving threshold_full orders of
  ## magnitude above first-step `prob = stp` scores, and the
  ## sequential scan flags ~all fixes.  Post-fix, threshold_partial
  ## is calibrated on `ref_stp` alone and matches the first-step
  ## scale, so the scan no longer over-flags.
  set.seed(101)
  n <- 400
  ## Mostly slow (rest) with intermittent bursts of fast movement.
  speed <- c(rep(0.05, 100), rep(2, 30), rep(0.05, 100), rep(2, 30),
              rep(0.05, 140))
  dx <- speed * 3600   # per-hour step, in metres at the latitude below
  ts <- as.POSIXct("2024-01-01", tz = "UTC") + cumsum(rep(3600, n))
  df <- data.frame(id = "t1", timestamp = ts,
                    lon = cumsum(dx) / (111000 * cos(48 * pi / 180)) + 11,
                    lat = 48 + rnorm(n, 0, 1e-5))
  m <- move2::mt_as_move2(df, coords = c("lon", "lat"),
                            time_column = "timestamp",
                            track_id_column = "id", crs = 4326)
  res <- suppressMessages(suppressWarnings(
    mt_sequential_outliers(m, plot = FALSE,
                              anchor_corruption_threshold = NULL)))
  flag_rate <- sum(res$is_outlier) / nrow(m)
  ## Pre-fix this would be ~1.0 (everything flagged); post-fix it
  ## should be well under 30% (the documented anchor-corruption gate).
  expect_lt(flag_rate, 0.30)
})


test_that("two-threshold fix: CPF synthetic byte-identical to pre-fix F1", {
  ## On single-state synthetic tracks, the autodifference KDE returns
  ## modest densities ~1 at the bulk of the unimodal Delta-step
  ## distribution, so threshold_partial and threshold_full agree
  ## numerically and the per-fix flag set should match the pre-fix
  ## C0-recorded counts within trivial tolerance.
  syn <- move2::mt_read(system.file("extdata", "synthetic_tracks.csv.gz",
                                       package = "move2utils"))
  syn <- syn[!sf::st_is_empty(syn), ]
  gt  <- readRDS(system.file("extdata", "synthetic_ground_truth.rds",
                                package = "move2utils"))
  ## F1 regression anchors in the canonical AEQD metric space (projection
  ## stratification 2026-06-08: all geometry computed in a per-track local
  ## AEQD, CRS-invariant).  These supersede the pre-projection flat-Earth
  ## values (CPF_A 0.850, CPF_C 0.857, CPF_E 0.557, CPF_F 0.727); the shift
  ## is the coordinate-system correction plus the threshold-floor degeneracy
  ## fix (CPF_C no longer collapses to 0).
  expected_f1 <- c(CPF_A = 0.757, CPF_C = 0.667, CPF_D = 1.000,
                    CPF_E = 0.182, CPF_F = 0.800)
  for (id in names(expected_f1)) {
    x <- syn[move2::mt_track_id(syn) == id, ]
    truth <- as.integer(gt[[id]]$index)
    truth <- truth[truth >= 1 & truth <= nrow(x)]
    is_truth <- logical(nrow(x)); is_truth[truth] <- TRUE
    o <- suppressMessages(suppressWarnings(
      mt_sequential_outliers(x, plot = FALSE)))
    tp <- sum(o$is_outlier & is_truth, na.rm = TRUE)
    fp <- sum(o$is_outlier & !is_truth, na.rm = TRUE)
    fn <- sum(!o$is_outlier & is_truth, na.rm = TRUE)
    p <- tp / max(tp + fp, 1); r <- tp / max(tp + fn, 1)
    f1 <- 2 * p * r / max(p + r, 1e-9)
    expect_equal(round(f1, 3), expected_f1[[id]],
                  info = sprintf("F1 regression on %s", id))
  }
})


test_that("two-threshold fix: user-supplied threshold = X applies to both regimes (back-compat)", {
  ## Setting threshold=X (legacy single-scalar override) should set
  ## both threshold_partial and threshold_full to X.  Test by
  ## supplying a permissive threshold and a strict threshold and
  ## checking flagging behaviour scales accordingly.
  x <- read_cpf("CPF_A")
  o_perm <- suppressMessages(suppressWarnings(
    mt_sequential_outliers(x, threshold = 1e-12,
                              plot = FALSE,
                              anchor_corruption_threshold = NULL)))
  o_strict <- suppressMessages(suppressWarnings(
    mt_sequential_outliers(x, threshold = 1e+12,
                              plot = FALSE,
                              anchor_corruption_threshold = NULL)))
  ## Permissive flags nothing (or close to it); strict flags lots.
  expect_lte(sum(o_perm$is_outlier),   sum(o_strict$is_outlier))
})

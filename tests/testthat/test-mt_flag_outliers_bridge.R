## Tests for mt_flag_outliers_bridge() and the shared threshold helpers.

suppressPackageStartupMessages({
  library(move2); library(sf)
})


# ---- fixtures ---------------------------------------------------------------

## Project a WGS84 move2 to a metric CRS using a local transverse-Mercator.
project_local <- function(x) {
  cc <- sf::st_coordinates(x)
  sf::st_transform(x, sprintf(
    "+proj=tmerc +lon_0=%f +lat_0=%f +ellps=WGS84 +units=m",
    mean(cc[, 1], na.rm = TRUE), mean(cc[, 2], na.rm = TRUE)))
}

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


# ---- basic input validation -------------------------------------------------

test_that("mt_flag_outliers_bridge rejects non-move2 input", {
  expect_error(mt_flag_outliers_bridge(data.frame(x = 1)),
                 "must be a move2 object")
})

test_that("mt_flag_outliers_bridge auto-projects lon/lat and returns in original CRS", {
  m   <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  ## the input is WGS84 (lon/lat); the function should auto-project,
  ## do the math, and return an object still in WGS84.
  expect_true(sf::st_is_longlat(m_C))
  res <- suppressMessages(
    mt_flag_outliers_bridge(m_C, threshold_type = "entropy", plot = FALSE))
  expect_identical(sf::st_crs(res), sf::st_crs(m_C))
  expect_equal(nrow(res), nrow(m_C))
  ## the flag columns should still be populated and in metres
  expect_true("bridge_residual" %in% names(res))
  expect_true(any(is.finite(res$bridge_residual)))
})

test_that("auto-projected and hand-projected inputs give identical flags", {
  m   <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ]
  m_C_proj <- project_local(m_C)
  res_auto <- suppressMessages(
    mt_flag_outliers_bridge(m_C, threshold_type = "entropy", plot = FALSE))
  res_hand <- suppressMessages(
    mt_flag_outliers_bridge(m_C_proj, threshold_type = "entropy", plot = FALSE))
  ## Same fixes should be flagged either way (tiny CRS choice difference is
  ## acceptable for η/break but the binary flag set should match).
  expect_equal(which(res_auto$is_outlier), which(res_hand$is_outlier))
})


# ---- recovery on CPF_C (golden test) ----------------------------------------

test_that("entropy threshold recovers all 4 CPF_C outliers with zero FPs", {
  m  <- read_synthetic()
  gt <- read_ground_truth()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])
  truth <- gt$CPF_C$index

  res <- suppressMessages(
    mt_flag_outliers_bridge(m_C, threshold_type = "entropy", plot = FALSE))
  flagged <- which(res$is_outlier)

  expect_setequal(flagged, truth)
  expect_equal(sum(res$is_outlier), length(truth))
})


# ---- precision on clean data ------------------------------------------------

test_that("entropy threshold flags no outliers on clean CPF_B", {
  m <- read_synthetic()
  m_B <- project_local(m[move2::mt_track_id(m) == "CPF_B", ])

  res <- suppressMessages(
    mt_flag_outliers_bridge(m_B, threshold_type = "entropy", plot = FALSE))
  expect_equal(sum(res$is_outlier), 0)
})


# ---- return object structure ------------------------------------------------

test_that("returned move2 carries expected diagnostic columns", {
  m <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])

  res <- suppressMessages(
    mt_flag_outliers_bridge(m_C, plot = FALSE))

  for (col in c("bridge_residual", "bridge_width", "bridge_eta",
                "bridge_percentile", "bridge_iteration",
                "is_outlier", "is_na_prob")) {
    expect_true(col %in% names(res),
                info = sprintf("column %s missing from result", col))
  }

  ## flag type sanity
  expect_type(res$is_outlier, "logical")
  expect_type(res$is_na_prob, "logical")

  ## first and last point have no bridge defined
  expect_true(is.na(res$bridge_eta[1]))
  expect_true(is.na(res$bridge_eta[nrow(res)]))
})


# ---- remove = TRUE drops flagged rows --------------------------------------

test_that("remove = TRUE returns object with flagged rows removed", {
  m  <- read_synthetic()
  gt <- read_ground_truth()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])

  res_keep <- suppressMessages(
    mt_flag_outliers_bridge(m_C, plot = FALSE, remove = FALSE))
  res_drop <- suppressMessages(
    mt_flag_outliers_bridge(m_C, plot = FALSE, remove = TRUE))

  expect_equal(nrow(res_drop),
                nrow(res_keep) - sum(res_keep$is_outlier))
  expect_false(any(res_drop$is_outlier))
})


# ---- multi-individual dispatch ---------------------------------------------

test_that("multi-individual input is processed per-track", {
  m <- project_local(read_synthetic())
  res <- suppressMessages(
    mt_flag_outliers_bridge(m, threshold_type = "entropy", plot = FALSE))

  expect_equal(nrow(res), nrow(m))
  expect_setequal(unique(as.character(move2::mt_track_id(res))),
                  c("CPF_A", "CPF_B", "CPF_C", "CPF_D", "CPF_E", "CPF_F"))
  ## bridge flags should appear in at least CPF_A and CPF_C
  per_id <- tapply(res$is_outlier, move2::mt_track_id(res), sum)
  expect_gt(per_id[["CPF_C"]], 0)
})


# ---- threshold_type argument validation -------------------------------------

test_that("threshold_type rejects unknown values", {
  m <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])
  expect_error(mt_flag_outliers_bridge(m_C, threshold_type = "bogus"),
                 "should be one of")
})


# ---- gap threshold finds more (with FPs) than entropy -----------------------

test_that("gap threshold is more sensitive than entropy on CPF_A", {
  m  <- read_synthetic()
  gt <- read_ground_truth()
  m_A <- project_local(m[move2::mt_track_id(m) == "CPF_A", ])

  res_e <- suppressMessages(
    mt_flag_outliers_bridge(m_A, threshold_type = "entropy", plot = FALSE))
  res_g <- suppressMessages(
    mt_flag_outliers_bridge(m_A, threshold_type = "gap", plot = FALSE))

  ## gap should flag at least as many as entropy (entropy is conservative)
  expect_gte(sum(res_g$is_outlier), sum(res_e$is_outlier))
})


# ---- short track degrades gracefully ---------------------------------------

test_that("too-short track returns unflagged with warning", {
  m  <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])
  m_small <- m_C[1:8, ]
  expect_warning(
    res <- mt_flag_outliers_bridge(m_small, plot = FALSE),
    "Too few")
  expect_equal(sum(res$is_outlier), 0)
})


# ---- internal helpers -------------------------------------------------------

test_that(".gap_threshold_lower returns consistent structure", {
  set.seed(1)
  ## bulk + spread-out outliers (mirrors real log(eta) distributions
  ## where outliers span orders of magnitude rather than clustering)
  outliers <- seq(-25, -10, length.out = 5) + rnorm(5, sd = 0.3)
  x <- c(rnorm(100), outliers)
  res <- move2utils:::.gap_threshold_lower(x)
  expect_named(res, c("is_outlier", "break_value", "method", "percentile"))
  expect_length(res$is_outlier, length(x))
  ## all injected outliers must be flagged
  expect_true(all(res$is_outlier[101:105]))
  ## the break must sit below the bulk median
  expect_lt(res$break_value, stats::median(x[1:100]))
})

test_that(".entropy_threshold_lower returns no outliers on unimodal data", {
  set.seed(1)
  x <- rnorm(500)  # unimodal, no valley
  res <- move2utils:::.entropy_threshold_lower(x)
  expect_equal(sum(res$is_outlier), 0)
  expect_identical(res$method, "none")
})

test_that(".entropy_threshold_lower finds valley on bimodal data", {
  set.seed(1)
  x <- c(rnorm(500, mean = 0), rnorm(20, mean = -5))
  res <- move2utils:::.entropy_threshold_lower(x)
  expect_gt(sum(res$is_outlier), 0)
  expect_identical(res$method, "valley")
})


# ---- directional method ------------------------------------------------------------

test_that("directional method recovers CPF_C ground truth with zero FPs", {
  m  <- read_synthetic()
  gt <- read_ground_truth()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])
  truth <- gt$CPF_C$index

  res <- suppressMessages(
    mt_flag_outliers_bridge(m_C, method = "directional", plot = FALSE))

  expect_setequal(which(res$is_outlier), truth)
})

test_that("directional method flags nothing on clean CPF_B", {
  m <- read_synthetic()
  m_B <- project_local(m[move2::mt_track_id(m) == "CPF_B", ])

  res <- suppressMessages(
    mt_flag_outliers_bridge(m_B, method = "directional", plot = FALSE))
  expect_equal(sum(res$is_outlier), 0)
})

test_that("all methods carry the full set of diagnostic eta columns", {
  m <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])

  for (meth in c("combined", "isotropic", "directional")) {
    res <- suppressMessages(
      mt_flag_outliers_bridge(m_C, method = meth, plot = FALSE))
    expect_true("bridge_eta"      %in% names(res))
    expect_true("bridge_eta_para" %in% names(res))
    expect_true("bridge_eta_perp" %in% names(res))
    ## magnitudes are non-negative
    expect_true(all(na.omit(res$bridge_eta_para) >= 0))
    expect_true(all(na.omit(res$bridge_eta_perp) >= 0))
    ## Pythagoras: eta^2 = eta_para^2 + eta_perp^2 (up to finite tol)
    ok <- is.finite(res$bridge_eta)
    expect_equal(res$bridge_eta[ok]^2,
                 (res$bridge_eta_para[ok]^2 + res$bridge_eta_perp[ok]^2),
                 tolerance = 1e-8)
  }
})

test_that("method argument rejects unknown values", {
  m <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])
  expect_error(mt_flag_outliers_bridge(m_C, method = "bogus"),
                 "should be one of")
})

test_that("directional preserves bridge_method attribute", {
  m <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])
  res <- suppressMessages(
    mt_flag_outliers_bridge(m_C, method = "directional", plot = FALSE))
  expect_identical(attr(res, "bridge_method"), "directional")
})

test_that("directional works on multi-individual input", {
  m <- project_local(read_synthetic())
  res <- suppressMessages(
    mt_flag_outliers_bridge(m, method = "directional", plot = FALSE))
  expect_equal(nrow(res), nrow(m))
  expect_true("bridge_eta_para" %in% names(res))
  expect_true("bridge_eta_perp" %in% names(res))
  ## multi-track should still flag CPF_C outliers
  per_id <- tapply(res$is_outlier, move2::mt_track_id(res), sum)
  expect_gt(per_id[["CPF_C"]], 0)
})

test_that("directional remove = TRUE drops flagged rows", {
  m <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])

  res_keep <- suppressMessages(
    mt_flag_outliers_bridge(m_C, method = "directional", plot = FALSE,
                             remove = FALSE))
  res_drop <- suppressMessages(
    mt_flag_outliers_bridge(m_C, method = "directional", plot = FALSE,
                             remove = TRUE))
  expect_equal(nrow(res_drop),
                nrow(res_keep) - sum(res_keep$is_outlier))
  expect_false(any(res_drop$is_outlier))
})


# ---- dirty-input hygiene ----------------------------------------------------
#
# Regression coverage for mixed hygiene failures that historically
# crashed .compute_bridge_residuals_dbgb() (NAs in subscript) and
# for the new hard-stop contract on duplicate / out-of-order times.

test_that("empty geometries are rejected with a classed error", {
  m <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])

  ## inject empty geometries
  ins_idx <- c(5L, 20L, 50L)
  geom <- sf::st_geometry(m_C)
  for (i in ins_idx) geom[[i]] <- sf::st_point()
  sf::st_geometry(m_C) <- geom

  ## As of 0.4.2 empty geometry is no longer silently skipped -- there is no
  ## outlier to identify in a fix with no location, so it must be removed
  ## upstream and the detectors reject it at entry.
  expect_error(
    suppressMessages(
      mt_flag_outliers_bridge(m_C, method = "combined", plot = FALSE)),
    class = "move2utils_input_has_empty_geometries")
})

test_that("duplicate timestamps hard-stop with mt_filter_unique guidance", {
  m <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])

  ## create a duplicate timestamp by cloning one row's time onto the next
  t <- move2::mt_time(m_C)
  t[10] <- t[9]
  attr(m_C, "time_column") <- "timestamp_dup"
  m_C$timestamp_dup <- t

  expect_error(
    suppressMessages(
      mt_flag_outliers_bridge(m_C, method = "combined", plot = FALSE)),
    "mt_filter_unique"
  )
})

test_that("out-of-order timestamps hard-stop with sort guidance", {
  m <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])

  ## swap two adjacent times to induce a negative Delta-t
  t <- move2::mt_time(m_C)
  tmp <- t[9]; t[9] <- t[10]; t[10] <- tmp
  attr(m_C, "time_column") <- "timestamp_swap"
  m_C$timestamp_swap <- t

  expect_error(
    suppressMessages(
      mt_flag_outliers_bridge(m_C, method = "combined", plot = FALSE)),
    "not time-sorted|negative Delta"
  )
})


# ---- location_error injection ---------------------------------------------------

test_that("location_error = NULL gives identical results to the default", {
  m  <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])
  res_default <- suppressMessages(
    mt_flag_outliers_bridge(m_C, plot = FALSE))
  res_null <- suppressMessages(
    mt_flag_outliers_bridge(m_C, location_error = NULL, plot = FALSE))
  expect_identical(res_default$is_outlier, res_null$is_outlier)
  expect_equal(res_default$bridge_eta, res_null$bridge_eta)
})

test_that("location_error always emits the bridge_obs_inflation column", {
  m   <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])
  res <- suppressMessages(
    mt_flag_outliers_bridge(m_C, plot = FALSE))
  expect_true("bridge_obs_inflation" %in% names(res))
  ## NULL location_error: inflation == 1 wherever it is defined
  inflation <- res$bridge_obs_inflation
  expect_true(all(is.na(inflation) | abs(inflation - 1) < 1e-10))
})

test_that("location_error scalar inflates the denominator and lowers eta", {
  m   <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])

  res_off <- suppressMessages(
    mt_flag_outliers_bridge(m_C, location_error = NULL, plot = FALSE))
  res_on  <- suppressMessages(
    mt_flag_outliers_bridge(m_C, location_error = 50, plot = FALSE))

  ## inflation > 1 on the active middle fixes
  infl <- res_on$bridge_obs_inflation
  expect_true(all(is.na(infl) | infl > 1 - 1e-9))
  expect_true(any(infl > 1 + 1e-6, na.rm = TRUE))

  ## eta on the bulk should drop (denominator inflates)
  ok <- is.finite(res_off$bridge_eta) & is.finite(res_on$bridge_eta)
  expect_lt(median(res_on$bridge_eta[ok]),
             median(res_off$bridge_eta[ok]))

  ## S_hat attribute attached when injection ran
  expect_true(!is.null(attr(res_on, "bridge_S_hat")))
  expect_true(is.finite(attr(res_on, "bridge_S_hat")))
})

test_that("location_error scalar preserves CPF_C truth recovery at 5 m sigma", {
  m   <- read_synthetic()
  gt  <- read_ground_truth()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])
  truth <- gt$CPF_C$index

  ## A small (5 m) per-fix obs-error is realistic for GPS and should
  ## not erode our truth-recovery on this benchmark track.
  res <- suppressMessages(
    mt_flag_outliers_bridge(m_C, location_error = 5,
                             threshold_type = "entropy", plot = FALSE))
  expect_setequal(which(res$is_outlier), truth)
})

test_that("location_error all-NA column is silently disabled", {
  m   <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])
  m_C$hacc <- NA_real_

  res <- suppressMessages(
    mt_flag_outliers_bridge(m_C, location_error = "hacc", plot = FALSE))
  ## injection disabled -> identical flags to NULL
  res_null <- suppressMessages(
    mt_flag_outliers_bridge(m_C, plot = FALSE))
  expect_identical(res$is_outlier, res_null$is_outlier)
  ## inflation should all be 1 (no injection)
  infl <- res$bridge_obs_inflation
  expect_true(all(is.na(infl) | abs(infl - 1) < 1e-10))
})

test_that("location_error = 'auto' returns NULL injection on synthetic data", {
  m   <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])
  ## synthetic_tracks has no eobs_hacc / argos_lc columns
  res <- suppressMessages(
    mt_flag_outliers_bridge(m_C, location_error = "auto", plot = FALSE))
  res_null <- suppressMessages(
    mt_flag_outliers_bridge(m_C, plot = FALSE))
  expect_identical(res$is_outlier, res_null$is_outlier)
})

test_that("location_error reads a user-supplied per-fix column", {
  m   <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])
  m_C$my_hacc <- rep(8, nrow(m_C))   # uniform 8 m sigma
  res <- suppressMessages(
    mt_flag_outliers_bridge(m_C, location_error = "my_hacc", plot = FALSE))
  ## inflation > 1 across the bulk
  expect_true(any(res$bridge_obs_inflation > 1 + 1e-6, na.rm = TRUE))
})

test_that("location_error rejects wrong-length numeric vector", {
  m   <- read_synthetic()
  m_C <- project_local(m[move2::mt_track_id(m) == "CPF_C", ])
  expect_error(
    suppressMessages(
      mt_flag_outliers_bridge(m_C,
                              location_error = c(5, 10),  # too short
                              plot = FALSE)),
    "must have length")
})

test_that("location_error suppresses sub-noise jitter in burst-sampled segment", {
  ## Construct a clean straight track, insert a burst-sampled segment
  ## (1 s gaps) where coordinates jitter by a few metres.  Without
  ## location_error injection, the bridge over-flags the jitter as
  ## outliers (small dt -> small bridge width -> large eta).  With a
  ## realistic location_error, those rate-extreme but absolute-tiny
  ## displacements should no longer be flagged.
  set.seed(2)
  n_main  <- 100
  n_burst <- 20
  ## main: 60 s sampling, 5 m/s along x with small BM noise so the
  ## residual scale S_hat has a defined positive value
  t_main <- seq(0, by = 60, length.out = n_main)
  x_main <- 5 * t_main + rnorm(n_main, 0, 6)
  y_main <- rnorm(n_main, 0, 6)
  ## burst inserted at t = 30000-30019 (1 s spacing); animal still
  ## moving at 5 m/s, but coordinates jittered with sd = 4 m
  t_burst <- 30000 + seq_len(n_burst)
  x_burst <- 5 * t_burst + rnorm(n_burst, 0, 4)
  y_burst <- rnorm(n_burst, 0, 4)

  t_all <- c(t_main[t_main < 30000], t_burst, t_main[t_main >= 30000])
  x_all <- c(x_main[t_main < 30000], x_burst, x_main[t_main >= 30000])
  y_all <- c(y_main[t_main < 30000], y_burst, y_main[t_main >= 30000])

  d <- data.frame(
    id = "T",
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") + t_all,
    x = x_all, y = y_all
  )
  m <- move2::mt_as_move2(d, coords = c("x", "y"),
                          time_column = "timestamp",
                          track_id_column = "id",
                          crs = "+proj=tmerc +lon_0=11 +lat_0=47 +ellps=WGS84 +units=m")

  res_off <- suppressMessages(
    mt_flag_outliers_bridge(m, location_error = NULL,
                             threshold_type = "gap", plot = FALSE))
  res_on  <- suppressMessages(
    mt_flag_outliers_bridge(m, location_error = 5,  # 5 m nominal GPS sigma
                             threshold_type = "gap", plot = FALSE))

  ## Mechanism check: where the bridge widths are narrowest (burst
  ## segment), the obs-error correction inflates the denominator most.
  ## Mean inflation in the burst block should be substantially higher
  ## than on the main track.
  burst_idx <- which(t_all >= 30000 & t_all < 30000 + n_burst)
  infl <- res_on$bridge_obs_inflation
  burst_mean <- mean(infl[burst_idx], na.rm = TRUE)
  main_mean  <- mean(infl[-burst_idx], na.rm = TRUE)
  expect_gt(burst_mean, main_mean * 2)

  ## Within the burst block, no more flags appear with the correction
  ## than without it (the absolute residual ~4 m is sub-noise once
  ## anchor sigma = 5 m enters the denominator).
  expect_lte(sum(res_on$is_outlier[burst_idx]),
              sum(res_off$is_outlier[burst_idx]))
})


# ---- .fn_core / wrapper-on-slice equivalence (cascade regression) ----------

## Regression test for the 2026-05-12 length-n vs length-k storage bug.
##
## Before the fix, `.bridge_fn_core` used length-n internal scratch with
## the active subset marked TRUE in a length-n mask.  This worked for
## the standalone wrapper path (where the wrapper called .fn_core with
## `active_idx_init = which(!bad_row)` -- almost always seq_len(n) on
## clean inputs) but silently diverged in the cascade hot path, where
## the cascade calls `.bridge_fn_core(cc_all, t_all_s, active_idx)`
## directly with `active_idx = which(!is_outlier)` shrinking each iter.
## With length-n internal scratch, cascade-already-flagged positions
## sat in the dedup arrays as FALSE -- breaking
## `.dedup_consecutive`'s "consecutive in active subset" semantics that
## the old wrapper-on-slice path produced (where those positions were
## absent from the sliced input entirely).  Two new flags separated by
## a cascade-already-flagged fix were kept-both under length-n
## storage; deduped to one under slice semantics.  Saline went from 32
## to 746 outliers as a result.
##
## The pre-existing test suite missed this because every test passed
## the bridge wrapper a contiguous track (no cascade-active-idx
## upstream).  This test calls `.bridge_fn_core` directly with a
## non-trivial active_idx and asserts it matches the wrapper-on-slice
## output (the historical cascade behaviour).

test_that(".bridge_fn_core matches wrapper-on-slice for non-trivial active_idx", {
  m   <- read_synthetic()
  ## Feed the canonical per-track AEQD (the projection the wrapper now
  ## canonicalises to) so the wrapper's projection is a no-op and the core
  ## (called on these same coords) matches it.  Under any other projection
  ## the wrapper would reproject and the two would legitimately differ.
  m_raw <- m[move2::mt_track_id(m) == "CPF_A", ]
  m_A <- sf::st_transform(
    m_raw, move2::mt_aeqd_crs(m_raw, center = "center", units = "m"))

  cc  <- sf::st_coordinates(m_A)
  t_s <- as.numeric(move2::mt_time(m_A), units = "secs")
  n   <- nrow(m_A)

  ## active_idx that skips every 13th fix -- many gaps between
  ## consecutive active positions, exercising the slice-semantics
  ## boundary where the bug previously manifested.
  active_idx <- setdiff(seq_len(n), seq(13L, n, by = 13L))

  res_core <- .bridge_fn_core(cc, t_s, active_idx,
                                method            = "combined",
                                threshold_type    = "entropy",
                                threshold         = 0.3,
                                residual_floor    = 0,
                                iterations        = 3L,
                                dedup_neighbours  = TRUE,
                                silent            = TRUE)

  res_wrap <- suppressMessages(
    mt_flag_outliers_bridge(m_A[active_idx, ],
                             plot = FALSE, silent = TRUE))

  ## is_outlier MUST be identical.  This is the primary regression
  ## assertion -- any divergence here means dedup semantics or active
  ## handling are out of sync between the cascade fast path and the
  ## wrapper slow path.
  expect_equal(res_core$is_outlier, res_wrap$is_outlier)

  ## bridge_residual should also match where both are non-NA.
  ## (Iter-by-iter shrinking can leave some positions never scored on
  ## one side; check the intersection of non-NA.)
  ok <- is.finite(res_core$bridge_residual) &
        is.finite(res_wrap$bridge_residual)
  expect_equal(res_core$bridge_residual[ok],
                res_wrap$bridge_residual[ok])
})


# ---- cascade multi-iter regression on cohort-like inputs -------------------

## ---- pool_by tests ---------------------------------------------------

## Helper: build a 2-track move2 sharing one animal id.
make_pair_bridge <- function(n_long = 200, n_short = 80, indv = "I1") {
  set.seed(7)
  build <- function(id, n, lon0, lat0) {
    dx <- rnorm(n, 0, 0.002); dy <- rnorm(n, 0, 0.002)
    data.frame(id = id,
               timestamp = as.POSIXct("2026-01-01", tz = "UTC") +
                            seq_len(n) * 3600,
               lon = lon0 + cumsum(dx), lat = lat0 + cumsum(dy))
  }
  d <- rbind(build("t1", n_long, 10, 50), build("t2", n_short, 12, 51))
  m <- move2::mt_as_move2(d, coords = c("lon", "lat"),
                            time_column = "timestamp",
                            track_id_column = "id", crs = 4326)
  td <- move2::mt_track_data(m); td$indv <- indv
  move2::mt_set_track_data(m, td)
}


test_that("bridge pool_by validates inputs", {
  m <- make_pair_bridge()
  expect_error(mt_flag_outliers_bridge(m, pool_by = 1, plot = FALSE),
                "length 1 \\(single column")
  expect_error(mt_flag_outliers_bridge(m, pool_by = c("a","b","c"),
                                          plot = FALSE),
                "Deeper hierarchies")
  expect_error(suppressMessages(
                  mt_flag_outliers_bridge(m, pool_by = "no_such_col",
                                            plot = FALSE, silent = TRUE)),
                "not in")
  expect_error(suppressMessages(
                  mt_flag_outliers_bridge(m, pool_by = c("indv","indv"),
                                            plot = FALSE, silent = TRUE)),
                "two \\*distinct\\* columns")
})


test_that("bridge pool_by = NULL is byte-identical to no pool_by", {
  m <- make_pair_bridge()
  for (mthd in c("combined", "isotropic", "directional")) {
    o1 <- suppressMessages(mt_flag_outliers_bridge(
            m, method = mthd, plot = FALSE, silent = TRUE))
    o2 <- suppressMessages(mt_flag_outliers_bridge(
            m, method = mthd, pool_by = NULL, plot = FALSE, silent = TRUE))
    expect_identical(o1$is_outlier, o2$is_outlier,
                      info = sprintf("method = %s", mthd))
    expect_identical(o1$bridge_eta, o2$bridge_eta,
                      info = sprintf("method = %s", mthd))
  }
})


test_that("bridge pool_by union is strictly additive (per-track ⊆ pool)", {
  ## Iterative bridge: per-track iteration uses track-local thresholds
  ## per-iter; pool fits one threshold on the converged lifted eta
  ## vector and re-applies.  The single-track-in-group byte-identity
  ## that holds for one-pass primitives (detour, speed_cap) does NOT
  ## hold for bridge in general -- the lifted eta is a mix of (final-
  ## iter eta for survivors) + (at-flag-iter eta for flagged fixes),
  ## so the pool break can differ from the per-track final-iter
  ## break.  What holds: pool flags are a superset of per-track
  ## flags (additive union contract).
  m <- make_pair_bridge(n_long = 300, n_short = 100)
  for (mthd in c("combined", "isotropic", "directional")) {
    o1 <- suppressMessages(mt_flag_outliers_bridge(
            m, method = mthd, plot = FALSE, silent = TRUE))
    o2 <- suppressMessages(mt_flag_outliers_bridge(
            m, method = mthd, pool_by = "indv",
            plot = FALSE, silent = TRUE))
    expect_true(all(which(o1$is_outlier) %in% which(o2$is_outlier)),
                 info = sprintf("method = %s: per-track NOT subset of pool",
                                 mthd))
  }
})


test_that("bridge pool union respects residual_floor gate", {
  ## A pool flag may only fire on a fix whose bridge_residual is
  ## strictly above the user-supplied residual_floor, mirroring the
  ## per-track gate.
  m <- make_pair_bridge(n_long = 200, n_short = 80)
  ## Set a very high residual_floor so nothing gets flagged.
  o <- suppressMessages(mt_flag_outliers_bridge(
          m, residual_floor = 1e9, pool_by = "indv",
          plot = FALSE, silent = TRUE))
  expect_equal(sum(o$is_outlier), 0L)
})


## Companion: directly assert mt_clean_track gives stable is_outlier
## counts on CPF tracks.  Locks in the post-refactor counts so any
## future refactor that silently shifts cleaning behaviour fails fast.

test_that("mt_clean_track is_outlier counts are stable on CPF synthetic", {
  m <- read_synthetic()
  ## Regression anchors for BOTH the canonical class_aware conjunction
  ## rule and the default evidence_corroborated rule (default flipped to
  ## the latter 2026-06-07; see DESIGN_evidence_accumulation.md sec 9c).
  ## counts in the canonical AEQD metric space (projection stratification
  ## 2026-06-08).  CPF_A ec was 24 under the test's project_local TMERC
  ## projection; canonicalisation now computes it in AEQD (= the longlat
  ## value, 25) -- the old 24 was the CRS-dependent artifact this fix removes.
  expected_ca <- c(CPF_A = 25L, CPF_B = 0L, CPF_C = 4L)  # conjunction + neighbours
  expected_ec <- c(CPF_A = 25L, CPF_B = 0L, CPF_C = 4L)  # evidence_corroborated (default)
  for (id in c("CPF_A", "CPF_B", "CPF_C")) {
    m_id <- project_local(m[move2::mt_track_id(m) == id, ])
    ca <- suppressMessages(suppressWarnings(
      mt_clean_track(m_id, plot = FALSE, remove = FALSE,
                     consensus = "class_aware")))
    ec <- suppressMessages(suppressWarnings(
      mt_clean_track(m_id, plot = FALSE, remove = FALSE)))  # default
    expect_equal(sum(ca$is_outlier), unname(expected_ca[id]),
                  info = sprintf("class_aware count on %s", id))
    expect_equal(sum(ec$is_outlier), unname(expected_ec[id]),
                  info = sprintf("evidence_corroborated (default) count on %s", id))
  }
})


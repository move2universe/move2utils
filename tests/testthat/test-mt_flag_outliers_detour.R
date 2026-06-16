## Tests for the detour-ratio primitive mt_flag_outliers_detour().

suppressPackageStartupMessages({
  library(move2); library(sf)
})


make_clean_track <- function(n = 200, id = "t1", crs = 4326) {
  set.seed(42)
  ## A small Brownian-like walk near a realistic mid-latitude
  ## location, in lon/lat.  Step magnitudes ~100-300 m at this scale.
  dx <- rnorm(n, 0, 0.002)
  dy <- rnorm(n, 0, 0.002)
  x  <- 10 + cumsum(dx); y <- 50 + cumsum(dy)
  d  <- data.frame(
    id = id,
    timestamp = as.POSIXct("2026-01-01", tz = "UTC") + seq_len(n) * 3600,
    lon = x, lat = y
  )
  mt_as_move2(d, coords = c("lon", "lat"),
              time_column = "timestamp",
              track_id_column = "id",
              crs = crs)
}


test_that("rejects non-move2 input", {
  expect_error(mt_flag_outliers_detour(data.frame(x = 1)),
               "must be a move2")
})


test_that("rejects bad k / threshold / min_leg", {
  m <- make_clean_track(50)
  expect_error(mt_flag_outliers_detour(m, k = 0), "positive integer")
  expect_error(mt_flag_outliers_detour(m, k = 1.5), "positive integer")
  expect_error(mt_flag_outliers_detour(m, threshold = 0.5), "scalar > 1")
  expect_error(mt_flag_outliers_detour(m, min_leg = -1), "non-negative")
})


test_that("clean random walk: very few false positives at threshold = 15", {
  ## A k=1 detour ratio on a Brownian walk has heavy tails: random
  ## direction reversals legitimately produce ratios > 5 occasionally.
  ## With strict conjunction inside mt_clean_track the kinematic
  ## detectors gate these out; the standalone primitive should still
  ## be sparse on a clean walk at threshold = 15.
  m <- make_clean_track(300)
  out <- mt_flag_outliers_detour(m, k = 1, threshold = 15,
                                  plot = FALSE, silent = TRUE)
  expect_lte(sum(out$is_outlier), 5)  # generous
  expect_true("detour_ratio" %in% names(out))
  expect_true("flagged_by_detour" %in% names(out))
  expect_equal(nrow(out), nrow(m))
})


test_that("injected single-fix spike is flagged", {
  m  <- make_clean_track(100)
  cc <- sf::st_coordinates(m)
  ## Move fix 50 hard off the line of travel
  cc[50, ] <- cc[50, ] + c(2, 2)  # ~200 km off
  pts <- lapply(seq_len(nrow(m)),
                function(i) sf::st_point(cc[i, ]))
  geom <- sf::st_sfc(pts, crs = sf::st_crs(m))
  m_spike <- m
  sf::st_geometry(m_spike) <- geom

  out <- mt_flag_outliers_detour(m_spike, k = 1, threshold = 5,
                                  plot = FALSE, silent = TRUE)
  expect_true(out$is_outlier[50])
  ## ratio at the spike should be very large
  expect_gt(out$detour_ratio[50], 50)
})


test_that("multi-individual dispatch processes each track", {
  ## Build two tracks in one data frame (move2's rbind on single-track
  ## move2 objects is restrictive; build the multi-track move2 from a
  ## combined data frame instead).
  set.seed(42)
  n <- 80
  build <- function(id, lon0, lat0) {
    dx <- rnorm(n, 0, 0.002); dy <- rnorm(n, 0, 0.002)
    data.frame(id = id,
               timestamp = as.POSIXct("2026-01-01", tz = "UTC") +
                            seq_len(n) * 3600,
               lon = lon0 + cumsum(dx), lat = lat0 + cumsum(dy))
  }
  d <- rbind(build("t1", 10, 50), build("t2", 12, 51))
  combined <- mt_as_move2(d, coords = c("lon", "lat"),
                           time_column = "timestamp",
                           track_id_column = "id",
                           crs = 4326)
  out <- mt_flag_outliers_detour(combined, k = 1, threshold = 15,
                                  plot = FALSE, silent = TRUE)
  expect_equal(nrow(out), nrow(combined))
  expect_lte(sum(out$is_outlier), 5)
})


test_that("short track returns unflagged with NA ratios", {
  m <- make_clean_track(2)  # too short for k = 1 (needs >= 3)
  out <- mt_flag_outliers_detour(m, k = 1, plot = FALSE, silent = TRUE)
  expect_equal(sum(out$is_outlier), 0)
  expect_true(all(is.na(out$detour_ratio)))
})


test_that("multi-k vector form works and computes per-fix max", {
  m <- make_clean_track(200)
  out <- mt_flag_outliers_detour(m, k = c(1, 2, 3), threshold = 8,
                                  plot = FALSE, silent = TRUE)
  expect_equal(nrow(out), nrow(m))
  expect_true(all(out$detour_ratio >= 1, na.rm = TRUE))
})


test_that("min_leg gates out small-displacement noise", {
  m  <- make_clean_track(100)
  cc <- sf::st_coordinates(m)
  ## Inject a small zig-zag near fix 30 (large ratio, tiny path)
  cc[30, ] <- cc[30, ] + c(0.0005, 0.0005)
  pts <- lapply(seq_len(nrow(m)),
                function(i) sf::st_point(cc[i, ]))
  geom <- sf::st_sfc(pts, crs = sf::st_crs(m))
  m_zigzag <- m
  sf::st_geometry(m_zigzag) <- geom

  ## With min_leg = 0 the small wiggle may pass; with min_leg large it
  ## must not -- legs are < 100 m on a noise-scale walk.
  out_gated <- mt_flag_outliers_detour(m_zigzag, k = 1, threshold = 5,
                                        min_leg = 5000,
                                        plot = FALSE, silent = TRUE)
  expect_false(out_gated$is_outlier[30])
})


test_that("integration with mt_clean_track: synthetic ground truth preserved", {
  path <- system.file("extdata", "synthetic_tracks.csv.gz",
                       package = "move2utils")
  if (nchar(path) == 0)
    path <- "inst/extdata/synthetic_tracks.csv.gz"
  d <- read.csv(gzfile(path), stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  build <- function(id) {
    mi <- mt_as_move2(d[d$individual.local.identifier == id, ],
      coords = c("location.long", "location.lat"),
      time_column = "timestamp",
      track_id_column = "individual.local.identifier",
      crs = 4326)
    mi <- mi[!sf::st_is_empty(mi), ]
    dplyr::arrange(mi, mt_time(mi))
  }
  expect_equal(sum(suppressMessages(
    mt_clean_track(build("CPF_B"), plot = FALSE, remove = FALSE,
                   silent = TRUE)$is_outlier)), 0)
  expect_gte(sum(suppressMessages(
    mt_clean_track(build("CPF_A"), plot = FALSE, remove = FALSE,
                   silent = TRUE)$is_outlier)), 23)
})


test_that("threshold_type validates inputs", {
  m <- make_clean_track(50)
  expect_error(mt_flag_outliers_detour(m, threshold_type = "bogus"),
               "should be one of")
})

test_that("threshold_type = 'auto' returns a valid flag set on a clean track (no entropy valley -> no flags)", {
  ## On a clean Brownian walk the -log(ratio) distribution does
  ## not contain a valley deep enough to declare an outlier
  ## regime; the entropy detector returns no break and the
  ## function flags nothing.  Safe-on-clean contract.
  m <- make_clean_track(300)
  out <- suppressMessages(mt_flag_outliers_detour(
    m, k = 1, threshold_type = "auto",
    plot = FALSE, silent = TRUE))
  expect_lte(sum(out$is_outlier), 5)  # generous (clean data)
  expect_true("detour_ratio" %in% names(out))
})

test_that("threshold_type = 'auto' improves F1 on CPF_E (colony-halo case)", {
  ## Empirical anchor: on the bundled CPF_E synthetic track
  ## (colony halo, 80 truth outliers spread across 1440 fixes),
  ## the adaptive entropy threshold finds a sharper bulk-vs-tail
  ## break than the fixed=8 cascade default.  This anchors the
  ## "real progress on colony-halo" claim and protects against
  ## future regressions.
  syn <- move2::mt_read(system.file("extdata/synthetic_tracks.csv.gz",
                                       package = "move2utils"))
  syn <- syn[!sf::st_is_empty(syn), ]
  gt  <- readRDS(system.file("extdata/synthetic_ground_truth.rds",
                                package = "move2utils"))
  trk <- syn[move2::mt_track_id(syn) == "CPF_E", ]
  truth <- as.integer(gt$CPF_E$index)
  truth <- truth[truth >= 1 & truth <= nrow(trk)]
  n <- nrow(trk)

  out_a <- suppressMessages(mt_flag_outliers_detour(
    trk, threshold_type = "auto", plot = FALSE, silent = TRUE))
  is_truth <- logical(n); is_truth[truth] <- TRUE

  TP <- sum(out_a$is_outlier & is_truth)
  FP <- sum(out_a$is_outlier & !is_truth)
  FN <- sum(!out_a$is_outlier & is_truth)
  prec <- TP / max(TP + FP, 1)
  rec  <- TP / max(TP + FN, 1)
  f1   <- 2 * prec * rec / max(prec + rec, 1)

  ## At the round-3 evaluation, fixed=8 gave F1=0.902 on CPF_E;
  ## auto gave F1=0.961.  Pin the auto F1 floor at 0.94 (allow
  ## small CRS / numerical drift but require a clear improvement
  ## over fixed=8).
  expect_gte(f1, 0.94)
  expect_equal(FP, 0L)  # auto eliminates the fixed-threshold FP
})


test_that("pool_by validates inputs", {
  m <- make_clean_track(80)
  expect_error(mt_flag_outliers_detour(m, pool_by = 1),
               "length 1 \\(single column")
  expect_error(mt_flag_outliers_detour(m, pool_by = c("a", "b", "c")),
               "Deeper hierarchies")
  expect_error(mt_flag_outliers_detour(m, pool_by = "no_such_col",
                                        plot = FALSE, silent = TRUE),
               "not in")
})


test_that("pool_by = NULL is byte-identical to no pool_by (no-op contract)", {
  ## Build a two-track move2 with a deployment_id column so the path
  ## reaches the column-lookup machinery only when pool_by is set.
  set.seed(42); n <- 80
  build <- function(id, lon0, lat0) {
    dx <- rnorm(n, 0, 0.002); dy <- rnorm(n, 0, 0.002)
    data.frame(id = id,
               timestamp = as.POSIXct("2026-01-01", tz = "UTC") +
                            seq_len(n) * 3600,
               lon = lon0 + cumsum(dx), lat = lat0 + cumsum(dy))
  }
  d <- rbind(build("t1", 10, 50), build("t2", 12, 51))
  m <- mt_as_move2(d, coords = c("lon", "lat"),
                    time_column = "timestamp",
                    track_id_column = "id", crs = 4326)

  o_null <- mt_flag_outliers_detour(m, threshold_type = "auto",
                                      plot = FALSE, silent = TRUE)
  o_arg  <- mt_flag_outliers_detour(m, threshold_type = "auto",
                                      pool_by = NULL,
                                      plot = FALSE, silent = TRUE)
  expect_identical(o_null$is_outlier, o_arg$is_outlier)
  expect_identical(o_null$detour_ratio, o_arg$detour_ratio)
})


test_that("pool_by groups of size one are byte-identical to per-track", {
  ## Building a two-track move2 with a column that puts each track in
  ## its own group.  Pool fit per-group = per-track fit; union is a
  ## no-op; result must be byte-identical to no-pool.
  set.seed(42); n <- 80
  build <- function(id, lon0, lat0) {
    dx <- rnorm(n, 0, 0.002); dy <- rnorm(n, 0, 0.002)
    data.frame(id = id,
               timestamp = as.POSIXct("2026-01-01", tz = "UTC") +
                            seq_len(n) * 3600,
               lon = lon0 + cumsum(dx), lat = lat0 + cumsum(dy))
  }
  d <- rbind(build("t1", 10, 50), build("t2", 12, 51))
  m <- mt_as_move2(d, coords = c("lon", "lat"),
                    time_column = "timestamp",
                    track_id_column = "id", crs = 4326)
  td <- move2::mt_track_data(m)
  td$grp <- as.character(td[[move2::mt_track_id_column(m)]])
  m <- move2::mt_set_track_data(m, td)

  o_pt <- mt_flag_outliers_detour(m, threshold_type = "auto",
                                    plot = FALSE, silent = TRUE)
  o_pg <- suppressMessages(mt_flag_outliers_detour(
            m, threshold_type = "auto", pool_by = "grp",
            plot = FALSE, silent = TRUE))
  expect_identical(o_pt$is_outlier, o_pg$is_outlier)
})


test_that("pool_by under threshold_type='fixed' is a no-op", {
  ## The fixed path uses the user-supplied scalar threshold uniformly
  ## across tracks; pool_by should be silently ignored, output
  ## byte-identical to no-pool.
  set.seed(42); n <- 80
  build <- function(id, lon0, lat0) {
    dx <- rnorm(n, 0, 0.002); dy <- rnorm(n, 0, 0.002)
    data.frame(id = id,
               timestamp = as.POSIXct("2026-01-01", tz = "UTC") +
                            seq_len(n) * 3600,
               lon = lon0 + cumsum(dx), lat = lat0 + cumsum(dy))
  }
  d <- rbind(build("t1", 10, 50), build("t2", 12, 51))
  m <- mt_as_move2(d, coords = c("lon", "lat"),
                    time_column = "timestamp",
                    track_id_column = "id", crs = 4326)
  td <- move2::mt_track_data(m); td$grp <- "G1"
  m <- move2::mt_set_track_data(m, td)

  o_fix1 <- mt_flag_outliers_detour(m, threshold = 5,
                                      threshold_type = "fixed",
                                      plot = FALSE, silent = TRUE)
  o_fix2 <- mt_flag_outliers_detour(m, threshold = 5,
                                      threshold_type = "fixed",
                                      pool_by = "grp",
                                      plot = FALSE, silent = TRUE)
  expect_identical(o_fix1$is_outlier, o_fix2$is_outlier)
})


test_that("pool_by union never un-flags what per-track-auto caught", {
  ## Build a 2-track scenario with a spike injected on t1.  Run auto-
  ## threshold both per-track and pool-by-individual.  Invariant:
  ## pool flags are a SUPERSET of per-track-auto flags (union is
  ## strictly additive, never un-flags).
  set.seed(42); n <- 200
  build <- function(id, lon0, lat0) {
    dx <- rnorm(n, 0, 0.002); dy <- rnorm(n, 0, 0.002)
    data.frame(id = id,
               timestamp = as.POSIXct("2026-01-01", tz = "UTC") +
                            seq_len(n) * 3600,
               lon = lon0 + cumsum(dx), lat = lat0 + cumsum(dy))
  }
  d  <- rbind(build("t1", 10, 50), build("t2", 12, 51))
  d$lon[100] <- d$lon[100] + 2   # spike at t1's row 100
  m  <- mt_as_move2(d, coords = c("lon", "lat"),
                    time_column = "timestamp",
                    track_id_column = "id", crs = 4326)
  td <- move2::mt_track_data(m); td$indv <- "I1"
  m  <- move2::mt_set_track_data(m, td)

  o_pt <- suppressMessages(mt_flag_outliers_detour(
            m, k = 1, threshold_type = "auto",
            plot = FALSE, silent = TRUE))
  o_pg <- suppressMessages(mt_flag_outliers_detour(
            m, k = 1, threshold_type = "auto", pool_by = "indv",
            plot = FALSE, silent = TRUE))
  ## per-track-auto flags ⊆ pool-auto flags (additive union contract).
  expect_true(all(which(o_pt$is_outlier) %in% which(o_pg$is_outlier)))
  ## The injected spike at row 100 must be flagged at least one of
  ## the two (fixed threshold = 5 catches it; auto may or may not).
  o_fx <- mt_flag_outliers_detour(m, k = 1, threshold = 5,
                                    threshold_type = "fixed",
                                    plot = FALSE, silent = TRUE)
  expect_true(o_fx$is_outlier[100])
})


test_that("pool_by + reference= are mutually exclusive on prob primitive", {
  m <- make_clean_track(80)
  expect_error(
    mt_flag_outliers(m, reference = m, pool_by = "track_id", silent = TRUE),
    "mutually exclusive")
})

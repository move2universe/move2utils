test_that("mt_peel_speed validates inputs", {
  expect_error(mt_peel_speed(NULL, v_max = 30), "move2")
  m <- readRDS(system.file("extdata", "synthetic_ground_truth.rds",
                             package = "move2utils"))[[1]]
  d <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")),
                 stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  x <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier", crs = 4326)
  x <- x[!sf::st_is_empty(x), ]
  xc <- x[move2::mt_track_id(x) == "CPF_C", ]

  expect_error(mt_peel_speed(xc, v_max = NA), "positive scalar")
  expect_error(mt_peel_speed(xc, v_max = -1), "positive scalar")
  expect_error(mt_peel_speed(xc, v_max = c(1, 2)), "positive scalar")
  expect_error(mt_peel_speed(xc, v_max = 30, max_iter = 0),
                "positive integer")
})

test_that("mt_peel_speed returns required columns and attributes", {
  d <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")),
                 stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  x <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier", crs = 4326)
  x <- x[!sf::st_is_empty(x), ]
  xc <- x[move2::mt_track_id(x) == "CPF_C", ]
  res <- suppressMessages(mt_peel_speed(xc, v_max = 30))
  expect_true(all(c("is_outlier", "peel_iteration", "step_speed") %in%
                     names(res)))
  expect_equal(attr(res, "v_max_used"), 30)
  expect_true(is.logical(attr(res, "converged")))
  expect_true(is.integer(attr(res, "n_peel_iterations")))
})

test_that("mt_peel_speed is idempotent on already-clean data", {
  ## CPF_C max step speed is ~0.6 m/s; a cap of 30 m/s peels nothing.
  d <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")),
                 stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  x <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier", crs = 4326)
  x <- x[!sf::st_is_empty(x), ]
  xc <- x[move2::mt_track_id(x) == "CPF_C", ]
  res1 <- suppressMessages(mt_peel_speed(xc, v_max = 30))
  expect_equal(sum(res1$is_outlier), 0L)
  expect_true(attr(res1, "converged"))
  ## Idempotent: second pass on survivors peels nothing more
  res2 <- suppressMessages(mt_peel_speed(res1, v_max = 30))
  expect_equal(sum(res2$is_outlier), 0L)
})

test_that("mt_peel_speed peels obvious outliers iteratively", {
  ## Inject a cluster of 3 outliers ~1000 m apart at fake location, then
  ## peel at v_max = 10.  Expect all 3 peeled; survivors max_speed <= 10.
  d <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")),
                 stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  x <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier", crs = 4326)
  x <- x[!sf::st_is_empty(x), ]
  xc <- x[move2::mt_track_id(x) == "CPF_C", ]
  ## Inject: shift 3 consecutive fixes by 0.1 degrees (~10 km)
  bad <- 50:52
  cc <- sf::st_coordinates(xc)
  cc[bad, 1] <- cc[bad, 1] + 1     # ~100 km offset
  new_geom <- sf::st_as_sf(as.data.frame(cc), coords = c("X", "Y"),
                            crs = sf::st_crs(xc))
  sf::st_geometry(xc) <- sf::st_geometry(new_geom)

  res <- suppressMessages(mt_peel_speed(xc, v_max = 10))
  ## Expect the 3 cluster fixes and possibly their immediate
  ## neighbours to be peeled (iterative propagation).
  expect_true(sum(res$is_outlier) >= 3L)
  expect_true(all(res$is_outlier[bad]))
  ## Post-peel: no survivor should have step_speed exceeding v_max
  s <- res$step_speed[!res$is_outlier]
  expect_true(all(is.na(s) | s <= 10))
})

test_that("mt_clean_track with v_max triggers pre-peel and flags cluster", {
  ## Same injection; expect mt_clean_track to flag the cluster when
  ## v_max = 10, and NOT to over-flag (< 10 % of track).
  d <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")),
                 stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  x <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier", crs = 4326)
  x <- x[!sf::st_is_empty(x), ]
  xc <- x[move2::mt_track_id(x) == "CPF_C", ]
  bad <- 50:52
  cc <- sf::st_coordinates(xc)
  cc[bad, 1] <- cc[bad, 1] + 1     # ~100 km offset -- large enough
                                   # that iterative peel propagation
                                   # removes all three within v_max=10
  new_geom <- sf::st_as_sf(as.data.frame(cc), coords = c("X", "Y"),
                            crs = sf::st_crs(xc))
  sf::st_geometry(xc) <- sf::st_geometry(new_geom)

  res <- suppressMessages(mt_clean_track(xc, v_max = 10, plot = FALSE,
                                            remove = FALSE))
  expect_true(all(res$is_outlier[bad]))
  expect_lt(mean(res$is_outlier), 0.15)
})

test_that("mt_peel_speed validates aux_scores", {
  d <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")),
                 stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  x <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier", crs = 4326)
  x <- x[!sf::st_is_empty(x), ]
  xc <- x[move2::mt_track_id(x) == "CPF_C", ]

  expect_error(
    mt_peel_speed(xc, v_max = 10, aux_scores = "not numeric"),
    "numeric vector"
  )
  expect_error(
    mt_peel_speed(xc, v_max = 10, aux_scores = c(1, 2, 3)),
    sprintf("nrow\\(x\\) = %d", nrow(xc))
  )
})

test_that("mt_peel_speed asymmetric mode preserves clean neighbours of a 1-fix spike", {
  ## Use clean CPF_B (no native outliers, regular sampling) and inject
  ## a single spike fix.  The symmetric peel removes the spike + both
  ## neighbours (3 fixes); the asymmetric peel with a perfect
  ## auxiliary score (highest at the spike) removes only the spike.
  d <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")),
                 stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  x <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier", crs = 4326)
  x <- x[!sf::st_is_empty(x), ]
  xc <- x[move2::mt_track_id(x) == "CPF_B", ]

  bad <- 50L  # single spike fix
  cc <- sf::st_coordinates(xc)
  cc[bad, 1] <- cc[bad, 1] + 1     # ~100 km offset spike
  new_geom <- sf::st_as_sf(as.data.frame(cc), coords = c("X", "Y"),
                            crs = sf::st_crs(xc))
  sf::st_geometry(xc) <- sf::st_geometry(new_geom)

  ## Symmetric peel: removes the spike and both neighbours.
  res_sym <- suppressMessages(mt_peel_speed(xc, v_max = 10))
  expect_true(res_sym$is_outlier[bad])
  expect_gte(sum(res_sym$is_outlier), 3L)  # spike + 2 neighbours

  ## Asymmetric peel with perfect oracle score (1 at the spike, 0 elsewhere):
  ## removes only the spike.
  oracle <- rep(0, nrow(xc))
  oracle[bad] <- 1
  res_asym <- suppressMessages(
    mt_peel_speed(xc, v_max = 10, aux_scores = oracle))
  expect_true(res_asym$is_outlier[bad])
  expect_equal(sum(res_asym$is_outlier), 1L)
  expect_false(res_asym$is_outlier[bad - 1L])
  expect_false(res_asym$is_outlier[bad + 1L])
  ## Survivor speeds clean.
  s <- res_asym$step_speed[!res_asym$is_outlier]
  expect_true(all(is.na(s) | s <= 10))
})

test_that("mt_peel_speed asymmetric mode equals symmetric when no aux_scores", {
  ## NULL aux_scores reproduces the symmetric default exactly.
  d <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")),
                 stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  x <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier", crs = 4326)
  x <- x[!sf::st_is_empty(x), ]
  xc <- x[move2::mt_track_id(x) == "CPF_C", ]

  res_default <- suppressMessages(mt_peel_speed(xc, v_max = 10))
  res_null    <- suppressMessages(
    mt_peel_speed(xc, v_max = 10, aux_scores = NULL))
  expect_equal(res_default$is_outlier, res_null$is_outlier)
})

test_that("mt_peel_speed asymmetric mode falls back to symmetric on NA scores", {
  ## When aux_scores has NA on either endpoint of an offending edge,
  ## the function flags both endpoints (symmetric fallback) for that
  ## edge so the peel still converges.
  d <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")),
                 stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  x <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier", crs = 4326)
  x <- x[!sf::st_is_empty(x), ]
  xc <- x[move2::mt_track_id(x) == "CPF_B", ]
  bad <- 50L
  cc <- sf::st_coordinates(xc)
  cc[bad, 1] <- cc[bad, 1] + 1
  new_geom <- sf::st_as_sf(as.data.frame(cc), coords = c("X", "Y"),
                            crs = sf::st_crs(xc))
  sf::st_geometry(xc) <- sf::st_geometry(new_geom)

  ## All-NA scores: every offending edge falls back to symmetric.
  na_scores <- rep(NA_real_, nrow(xc))
  res_na <- suppressMessages(
    mt_peel_speed(xc, v_max = 10, aux_scores = na_scores))
  expect_true(res_na$is_outlier[bad])
  expect_gte(sum(res_na$is_outlier), 3L)
})

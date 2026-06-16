## Helper: build a deterministic move2 fixture with two tracks --
## (A) a 6-segment "corridor": net east-going but with back-and-forth
##     oscillation along the east axis. Pseudo-azimuth collapses 90 vs
##     270 onto the same value, so all six segments are directionally
##     consistent. The oscillation produces variable segment lengths,
##     ensuring adjacent segments fall within each other's half-segment
##     search radius (a perfectly-uniform straight track structurally
##     does not, in this algorithm).
## (B) a 6-segment random-direction "resting" cluster of comparable
##     fix count, well separated from A. Low speed, high directional
##     variance, no overlap with A.
## CRS = UTM31N so distances/speeds are exact in metres.
make_corridor_fixture <- function() {
  base_lon <- 600000  # easting, metres
  base_lat <- 4600000 # northing, metres
  ## (A) east-going corridor with east-axis oscillation: fixes at
  ##     0, 100, 50, 150, 100, 200, 150 m east of base.
  a_e <- base_lon + c(0, 100, 50, 150, 100, 200, 150)
  a_n <- rep(base_lat, 7L)
  ## (B) 7 fixes resting: tight random scatter within a 10 m circle,
  ##     2 km east of A, 100 s apart -> low speed and high variance.
  set.seed(42L)
  b_e <- base_lon + 5000 + stats::runif(7L, -5, 5)
  b_n <- base_lat        + stats::runif(7L, -5, 5)
  t0  <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  df  <- data.frame(
    track = rep(c("A", "B"), each = 7L),
    time  = rep(t0 + seq(0, 6) * 100, 2L),
    e     = c(a_e, b_e),
    n     = c(a_n, b_n)
  )
  sf_obj <- sf::st_as_sf(df, coords = c("e", "n"), crs = 32631)
  move2::mt_as_move2(sf_obj, time_column = "time",
                     track_id_column = "track")
}

test_that("mt_corridor rejects non-move2 input with a classed error", {
  expect_error(
    mt_corridor(data.frame(x = 1:3, y = 1:3)),
    class = "move2utils_input_not_move2"
  )
})

test_that("mt_corridor rejects input with fewer than 2 non-empty rows", {
  fx <- make_corridor_fixture()
  one <- fx[1L, ]
  expect_error(
    mt_corridor(one),
    class = "move2utils_mt_corridor_too_few_rows"
  )
})

test_that("mt_corridor rejects input with fewer than 2 segments", {
  ## one fix per track -> two POINTs, no LINESTRINGs.
  fx <- make_corridor_fixture()
  two <- fx[c(1L, 8L), ]  # one fix from each track
  expect_error(
    mt_corridor(two),
    class = "move2utils_mt_corridor_too_few_segments"
  )
})

test_that("mt_corridor accepts projected (non-longlat) input", {
  fx <- make_corridor_fixture()
  out <- suppressWarnings(mt_corridor(fx))
  expect_s3_class(out, "move2")
  expect_true(all(c("corridor", "corridor_speed", "corridor_azimuth",
                    "corridor_circvar", "corridor_n_neighbours")
                  %in% names(out)))
  expect_equal(nrow(out), nrow(fx))
})

test_that("mt_corridor accepts longlat input (AEQD reprojection path)", {
  fx <- make_corridor_fixture()
  fx_ll <- sf::st_transform(fx, 4326)
  out <- suppressWarnings(mt_corridor(fx_ll))
  expect_s3_class(out, "move2")
  expect_equal(sf::st_crs(out), sf::st_crs(fx_ll))
  expect_equal(nrow(out), nrow(fx_ll))
})

test_that("mt_corridor flags the straight track and not the resting cluster", {
  fx <- make_corridor_fixture()
  ## explicit thresholds make the assertion independent of the
  ## within-object quantile fallback.
  out <- mt_corridor(fx,
                     speed_threshold   = 0.5,    # m/s; A=1, B<<0.5
                     circvar_threshold = 0.05,   # A near 0, B near 1
                     min_segments      = 2L)
  ## track A: all interior segments should be flagged corridor.
  ## (the first/last segments may have asymmetric neighbour counts;
  ## interior segments are the strongest test.)
  a_rows <- which(move2::mt_track_id(out) == "A")
  a_interior <- a_rows[3:5]
  expect_true(all(out$corridor[a_interior] == "corridor"))
  ## track B: no segment should be flagged corridor.
  b_rows <- which(move2::mt_track_id(out) == "B")
  expect_true(all(out$corridor[b_rows] == "not corridor"))
})

test_that("mt_corridor per-segment columns match expected values on kept rows", {
  fx <- make_corridor_fixture()
  out <- suppressWarnings(mt_corridor(fx))
  expected_speed <- as.numeric(move2::mt_speed(fx, units = "m/s"))
  expect_equal(out$corridor_speed, expected_speed)
  ## Azimuth: corridor_azimuth is the planar bearing (atan2 on segment
  ## endpoints in the work CRS) in degrees, with NA at track ends.
  ## On the east-east-then-back fixture, segment 1 goes (+100, 0) east
  ## (bearing 90); segment 2 goes (-50, 0) west (bearing -90).
  expect_equal(out$corridor_azimuth[1L],  90, tolerance = 1e-6)
  expect_equal(out$corridor_azimuth[2L], -90, tolerance = 1e-6)
  expect_true(all(out$corridor_azimuth[1:6] %in% c(-90, 90)))
})

test_that("mt_corridor track-final rows carry NA for per-segment columns", {
  fx <- make_corridor_fixture()
  out <- suppressWarnings(mt_corridor(fx))
  tail_rows <- c(7L, 14L)  # last fix of each track in the fixture
  expect_true(all(is.na(out$corridor_speed[tail_rows])))
  expect_true(all(is.na(out$corridor_azimuth[tail_rows])))
})

test_that("mt_corridor does not pollute the neighbour search across tracks", {
  ## Place two tracks so the spatial-distance check would invite
  ## cross-track segments if mt_segments were not track-aware:
  ## a single fix of B is dropped into the middle of A's path.
  fx <- make_corridor_fixture()
  ## Build a control input where B's fixes overlap A's coordinates
  ## exactly; with proper track-awareness B is still treated as its
  ## own track and never produces a segment linking it to A.
  set.seed(0L)
  base_lon <- 600000; base_lat <- 4600000
  e <- c(base_lon + seq(0, 6) * 100, base_lon + seq(0, 6) * 100 + 50)
  n <- rep(base_lat, 14L)
  t0 <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  df <- data.frame(
    track = rep(c("A", "B"), each = 7L),
    time  = rep(t0 + seq(0, 6) * 100, 2L),
    e = e, n = n
  )
  sf_obj <- sf::st_as_sf(df, coords = c("e", "n"), crs = 32631)
  fx2 <- move2::mt_as_move2(sf_obj, time_column = "time",
                            track_id_column = "track")
  out <- suppressWarnings(mt_corridor(fx2))
  ## All segments are real LINESTRINGs (track-end rows carry NA);
  ## row 7 is A's last fix, row 14 is B's last fix.
  expect_true(is.na(out$corridor_speed[7L]))
  expect_true(is.na(out$corridor_speed[14L]))
  ## Track A's first segment (row 1) should have a non-NA speed
  ## computed only from A's own next fix.
  expect_false(is.na(out$corridor_speed[1L]))
  expect_equal(out$corridor_speed[1L],
               as.numeric(move2::mt_speed(fx2, units = "m/s"))[1L])
})

test_that("mt_corridor handles empty geometries with an inform + NA-pad", {
  fx <- make_corridor_fixture()
  ## Inject two empty POINT geometries (e.g. ACC-only rows).
  empty_pt <- sf::st_sfc(sf::st_point(), crs = sf::st_crs(fx))
  sf::st_geometry(fx)[c(3L, 10L)] <- empty_pt
  ## The notice is a classed message.
  expect_message(
    suppressWarnings(mt_corridor(fx)),
    class = "move2utils_mt_corridor_empty_dropped"
  )
  out <- suppressMessages(suppressWarnings(mt_corridor(fx)))
  expect_equal(nrow(out), nrow(fx))
  ## Empty rows are NA-padded on all 5 new columns.
  expect_true(is.na(out$corridor_speed[3L]))
  expect_true(is.na(out$corridor_circvar[10L]))
  ## Their corridor label defaults to "not corridor" (not NA).
  expect_equal(as.character(out$corridor[c(3L, 10L)]),
               c("not corridor", "not corridor"))
})

test_that("mt_corridor warns when speed_threshold / circvar_threshold default", {
  fx <- make_corridor_fixture()
  expect_warning(
    mt_corridor(fx, circvar_threshold = 0.05),
    class = "move2utils_mt_corridor_default_speed_threshold"
  )
  expect_warning(
    mt_corridor(fx, speed_threshold = 0.5),
    class = "move2utils_mt_corridor_default_circvar_threshold"
  )
})

test_that("mt_corridor warns on irregular sampling", {
  ## Build a track with deliberately ragged time intervals.
  base_lon <- 600000; base_lat <- 4600000
  e <- base_lon + seq(0, 6) * 100
  n <- rep(base_lat, 7L)
  t0 <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  ## intervals: 100, 100, 10000, 100, 100, 10000 seconds.
  times <- t0 + c(0, 100, 200, 10200, 10300, 10400, 20400)
  df <- data.frame(
    track = rep("A", 7L), time = times, e = e, n = n
  )
  sf_obj <- sf::st_as_sf(df, coords = c("e", "n"), crs = 32631)
  fx_irr <- move2::mt_as_move2(sf_obj, time_column = "time",
                               track_id_column = "track")
  expect_warning(
    mt_corridor(fx_irr,
                speed_threshold   = 0.5,
                circvar_threshold = 0.05),
    class = "move2utils_mt_corridor_irregular_sampling"
  )
})

test_that("mt_corridor verbose flag controls the per-call summary message", {
  fx <- make_corridor_fixture()
  ## verbose = FALSE (default) + explicit thresholds: no output.
  expect_silent(
    mt_corridor(fx, speed_threshold = 0.5, circvar_threshold = 0.05)
  )
  ## verbose = TRUE: message about corridor count.
  expect_message(
    mt_corridor(fx, speed_threshold = 0.5, circvar_threshold = 0.05,
                verbose = TRUE),
    "Found [0-9]+ corridor segments"
  )
})

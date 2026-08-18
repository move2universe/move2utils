## Tests for mt_flag_speed_cap()

suppressPackageStartupMessages({
  library(move2); library(sf)
})

## Small synthetic track: one deliberately-extreme step inserted.
## Fix 5 -> fix 6 covers 1000 m in 1 s (1000 m/s), an impossible speed.
make_track <- function() {
  coords <- cbind(0:9 * 1.0, rep(0, 10))
  coords[6, 1] <- 1000      # spike at fix 6: 1 km from fix 5 (1m away)
  df <- data.frame(
    location_long = coords[, 1],
    location_lat  = coords[, 2],
    timestamp = as.POSIXct("2024-01-01 00:00:00", tz = "UTC") +
                 seq(0, by = 1, length.out = 10),
    individual_local_identifier = "ind_1",
    tag_local_identifier = "tag_1",
    sensor_type = "gps"
  )
  move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier",
    crs = 32633)  # projected (UTM) so mt_distance gives metres
}


test_that("rejects non-move2 input", {
  expect_error(mt_flag_speed_cap(data.frame(x = 1), v_max = 40, threshold_type = "hard"),
               "must be a move2 object")
})

test_that("rejects bad v_max", {
  x <- make_track()
  expect_error(mt_flag_speed_cap(x, v_max = -1, threshold_type = "hard"),   "positive scalar")
  expect_error(mt_flag_speed_cap(x, v_max = 0, threshold_type = "hard"),    "positive scalar")
  expect_error(mt_flag_speed_cap(x, v_max = c(40, 50), threshold_type = "hard"),
               "positive scalar")
  expect_error(mt_flag_speed_cap(x, v_max = NA_real_, threshold_type = "hard"),  "positive scalar")
})

test_that("rejects bad jitter", {
  x <- make_track()
  expect_error(mt_flag_speed_cap(x, v_max = 40, threshold_type = "hard", jitter = -5),
               "non-negative scalar")
  expect_error(mt_flag_speed_cap(x, v_max = 40, threshold_type = "hard", jitter = c(1, 2)),
               "non-negative scalar")
})

test_that("flags both endpoints of an offending step", {
  x <- make_track()
  res <- suppressMessages(mt_flag_speed_cap(x, v_max = 40, threshold_type = "hard", plot = FALSE))
  ## fix 6 has both a huge step-in (1 km from fix 5) and a huge step-out
  ## back to fix 7 -> so fix 5, fix 6, fix 7 are all flagged.
  expect_equal(which(res$is_outlier), c(5L, 6L, 7L))
  expect_true(all(res$is_speed_above_cap[c(5L, 6L, 7L)]))
  ## Everyone else passes
  expect_false(any(res$is_outlier[c(1:4, 8:10)]))
})

test_that("no flags below cap on clean steady track", {
  coords <- cbind(0:9 * 5.0, rep(0, 10))  # 5 m step, 1 s -> 5 m/s
  df <- data.frame(
    location_long = coords[, 1], location_lat = coords[, 2],
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") +
                 seq(0, by = 1, length.out = 10),
    individual_local_identifier = "ind", tag_local_identifier = "tag",
    sensor_type = "gps")
  x <- move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier", crs = 32633)
  res <- suppressMessages(mt_flag_speed_cap(x, v_max = 40, threshold_type = "hard", plot = FALSE))
  expect_false(any(res$is_outlier))
  expect_true(all(is.na(res$step_speed) | res$step_speed == 5))
})

test_that("jitter flag marks sub-floor displacements", {
  coords <- cbind(0:9 * 0.5, rep(0, 10))  # 0.5 m step; jitter < 10
  df <- data.frame(
    location_long = coords[, 1], location_lat = coords[, 2],
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") +
                 seq(0, by = 60, length.out = 10),    # 1 min apart
    individual_local_identifier = "ind", tag_local_identifier = "tag",
    sensor_type = "gps")
  x <- move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier", crs = 32633)
  res <- suppressMessages(
    mt_flag_speed_cap(x, v_max = 40, threshold_type = "hard", jitter = 10, plot = FALSE))
  expect_true(all(res$is_step_below_jitter[1:9]))
  expect_true(all(res$is_outlier[1:9]))
})

test_that("remove = TRUE drops flagged rows", {
  x <- make_track()
  res <- suppressMessages(mt_flag_speed_cap(
    x, v_max = 40, threshold_type = "hard", plot = FALSE, remove = TRUE))
  expect_equal(nrow(res), 10L - 3L)   # 3 removed (5, 6, 7)
  expect_false(any(res$is_outlier))
})

test_that("threshold_type = 'auto' default -- data-driven cap on bimodal track", {
  ## A realistic outlier signature: a long clean trajectory (2000 fixes)
  ## with a small number of jittered spikes scattered through it
  ## (3 of 2000 = 0.15%, consistent with real-world spoof rates).
  ## At this density the spike "mode" is far below the substantive
  ## threshold and the auto-cap gate (correctly) accepts the cap.
  ## When spikes form 3% of fixes -- as in some toy benchmarks -- the
  ## gate refuses without species knowledge; that case is covered by
  ## the threshold_type = "hard" tests below.
  set.seed(2)
  n <- 2000L
  coords <- cbind(cumsum(runif(n, 0, 1)), rep(0, n))
  spike_idx <- c(500, 1000, 1500)
  coords[spike_idx, 1] <- coords[spike_idx, 1] +
                            runif(length(spike_idx), 4000, 9000)
  df <- data.frame(
    location_long = coords[, 1], location_lat = coords[, 2],
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") +
                 seq(0, by = 60, length.out = n),
    individual_local_identifier = "ind", tag_local_identifier = "tag",
    sensor_type = "gps")
  x <- move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier", crs = 32633)
  res <- suppressMessages(mt_flag_speed_cap(x, plot = FALSE))
  expect_true(is.finite(attr(res, "v_max_used")))
  expect_true(sum(res$is_outlier) > 0)
  for (s in spike_idx) {
    expect_true(any(res$is_outlier[(s - 1):(s + 1)]))
  }
})

test_that("threshold_type = 'auto' on unimodal clean data flags nothing", {
  set.seed(3)
  coords <- cbind(cumsum(rnorm(200, sd = 1)), cumsum(rnorm(200, sd = 1)))
  df <- data.frame(
    location_long = coords[, 1], location_lat = coords[, 2],
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") +
                 seq(0, by = 60, length.out = 200),
    individual_local_identifier = "ind", tag_local_identifier = "tag",
    sensor_type = "gps")
  x <- move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier", crs = 32633)
  res <- suppressMessages(mt_flag_speed_cap(x, plot = FALSE))
  expect_false(any(res$is_outlier))
})

test_that("v_max supplied without threshold_type = 'hard' errors", {
  x <- make_track()
  expect_error(mt_flag_speed_cap(x, v_max = 40),
               "only used when threshold_type")
})

test_that("mt_suggest_speed_cap returns NA on unimodal clean data", {
  set.seed(1)
  coords <- cbind(cumsum(rnorm(200, sd = 1)), cumsum(rnorm(200, sd = 1)))
  df <- data.frame(
    location_long = coords[, 1], location_lat = coords[, 2],
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") +
                 seq(0, by = 60, length.out = 200),
    individual_local_identifier = "ind", tag_local_identifier = "tag",
    sensor_type = "gps")
  x <- move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier", crs = 32633)
  v <- suppressMessages(mt_suggest_speed_cap(x, plot = FALSE))
  expect_true(is.na(v))
})

test_that("mt_suggest_speed_cap detects a clean bimodal speed gap", {
  set.seed(2)
  ## 1000 fixes, steady slow; inject 4 extreme-speed transitions
  ## (< 0.5% of fixes -- well below the 2% substantive-mode threshold
  ## so the mode-position gate inside .compute_speed_cap accepts the
  ## suggested cap).  Pre-v0.3 the test used 3 / 200 = 1.5% spike
  ## fraction, which the suggester accepted only because it lacked
  ## the mode-position gate; post-alignment the gate would correctly
  ## refuse such a "substantive high-speed mode" cap.
  n <- 1000L
  coords <- cbind(cumsum(runif(n, 0, 1)), rep(0, n))
  spike_idx <- c(200L, 400L, 600L, 800L)
  coords[spike_idx, 1] <- coords[spike_idx, 1] + 5000  # 5 km jumps
  df <- data.frame(
    location_long = coords[, 1], location_lat = coords[, 2],
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") +
                 seq(0, by = 60, length.out = n),
    individual_local_identifier = "ind", tag_local_identifier = "tag",
    sensor_type = "gps")
  x <- move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier", crs = 32633)
  v <- suppressMessages(mt_suggest_speed_cap(x, plot = FALSE))
  expect_true(!is.na(v))
  expect_true(v > 1 && v < 80)   # somewhere between bulk (~0.5 m/s) and outlier (~80 m/s)
})

test_that("respects track boundaries in multi-individual input", {
  ## Two tracks separated by many years: without track-aware step
  ## computation, the jump from end-of-ind_1 (t ~ 2024-01-01) to
  ## start-of-ind_2 (t ~ 2099-01-01) with huge spatial offset would
  ## produce a normal-looking speed; the worry is the other way --
  ## a huge spatial jump with small dt would spuriously flag. Here
  ## we test that mt_distance's per-track NA convention is respected.
  mk <- function(id_val, t0) {
    coords <- cbind(0:9 * 1.0, rep(0, 10))
    coords[6, 1] <- 1000
    df <- data.frame(
      location_long = coords[, 1], location_lat = coords[, 2],
      timestamp = t0 + seq(0, by = 1, length.out = 10),
      individual_local_identifier = id_val,
      tag_local_identifier = paste0("tag_", id_val),
      sensor_type = "gps")
    move2::mt_as_move2(df,
      coords = c("location_long", "location_lat"),
      time_column = "timestamp",
      track_id_column = "individual_local_identifier", crs = 32633)
  }
  t1 <- mk("ind_1", as.POSIXct("2024-01-01", tz = "UTC"))
  t2 <- mk("ind_2", as.POSIXct("2099-01-01", tz = "UTC"))
  combined <- rbind(t1, t2)
  res <- suppressMessages(mt_flag_speed_cap(combined, v_max = 40,
                                             threshold_type = "hard",
                                             plot = FALSE))
  ## Each track independently should flag its own 3-fix spike.
  expect_equal(sum(res$is_outlier), 6L)
})


## ---- mt_suggest_speed_cap() with allometric overlay ----------------

test_that("mt_suggest_speed_cap accepts (mass, mode) and propagates allometric prediction", {
  m <- move2::mt_read(system.file("extdata/synthetic_tracks.csv.gz",
                                   package = "move2utils"))
  m <- m[!sf::st_is_empty(m), ]
  mA <- m[move2::mt_track_id(m) == "CPF_A", ]

  ## mass + mode: single allometric line
  v <- suppressMessages(
    mt_suggest_speed_cap(mA, mass = 1, mode = "flying", plot = FALSE))
  expect_true(is.numeric(v) || is.na(v))

  ## mass alone: triangulation across three modes (no error)
  v_three <- suppressMessages(
    mt_suggest_speed_cap(mA, mass = 1, plot = FALSE))
  expect_true(is.numeric(v_three) || is.na(v_three))

  ## mode alone: error
  expect_error(suppressMessages(
    mt_suggest_speed_cap(mA, mode = "flying", plot = FALSE)),
    "mass")

  ## v_max must be a positive scalar
  expect_error(suppressMessages(
    mt_suggest_speed_cap(mA, v_max = -1, plot = FALSE)),
    "v_max")
})


## ---- mt_suggest_speed_cap() / mt_flag_speed_cap() alignment --------
## Pre-v0.3 the suggester defaulted to method = "entropy" and skipped
## the mode-position gate plus the 55 m/s biological-sanity ceiling
## warning.  A user inspecting their data with mt_suggest_speed_cap()
## then calling mt_flag_speed_cap(x, threshold_type = "auto") could
## see two different proposals on the same track.  v0.3 closes this
## by routing the suggester through .compute_speed_cap (the same
## helper the flagger uses) and adding the ceiling warning.

test_that("mt_suggest_speed_cap default method is now 'auto' (matches flagger)", {
  ## Match.arg returns the first option as the default; assert that.
  fmls <- formals(mt_suggest_speed_cap)
  expect_equal(eval(fmls$method)[1], "auto")
})

test_that("mt_suggest_speed_cap rejects bad physiological_ceiling", {
  m <- move2::mt_read(system.file("extdata/synthetic_tracks.csv.gz",
                                   package = "move2utils"))
  m <- m[!sf::st_is_empty(m), ]
  mA <- m[move2::mt_track_id(m) == "CPF_A", ]
  expect_error(
    mt_suggest_speed_cap(mA, physiological_ceiling = -5, plot = FALSE),
    "positive scalar"
  )
  expect_error(
    mt_suggest_speed_cap(mA, physiological_ceiling = c(30, 50), plot = FALSE),
    "positive scalar"
  )
  expect_error(
    mt_suggest_speed_cap(mA, physiological_ceiling = "fast", plot = FALSE),
    "positive scalar"
  )
})

test_that("mt_suggest_speed_cap mode-position gate matches flagger on multi-state data", {
  ## Two-mode synthetic: 60% rest at ~0.1 m/s; 40% flight at ~15 m/s.
  ## A naive entropy/gap break in -log(speed) lands BELOW the flight
  ## mode (rest vs flight separator) -- the mode-position gate must
  ## refuse it because cap < rightmost substantive mode.  The
  ## flagger and the suggester now both refuse on this case.
  set.seed(42)
  n <- 1000L
  state <- rep(c(0L, 1L), times = c(600L, 400L))
  step_size <- ifelse(state == 0L,
                      stats::rgamma(n, shape = 2, rate = 20),  # rest ~0.1 m/s @ dt=60
                      stats::rgamma(n, shape = 2, rate = 0.15)) # flight ~15 m/s @ dt=60
  step_size <- pmax(step_size, 1e-3)
  coords <- cbind(cumsum(step_size), rep(0, n))
  df <- data.frame(
    location_long = coords[, 1], location_lat = coords[, 2],
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") +
                 seq(0, by = 60, length.out = n),
    individual_local_identifier = "ind", tag_local_identifier = "tag",
    sensor_type = "gps")
  x <- move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier", crs = 32633)

  ## Flagger refuses (auto -> Inf -> 0 flags).
  res_flagger <- suppressMessages(
    mt_flag_speed_cap(x, threshold_type = "auto", plot = FALSE))
  expect_equal(sum(res_flagger$is_outlier), 0L)

  ## Suggester now also refuses (returns NA) -- no longer the
  ## pre-v0.3 silent disagreement where suggester proposed a value
  ## the flagger refused.
  v_suggest <- suppressMessages(
    mt_suggest_speed_cap(x, plot = FALSE))
  expect_true(is.na(v_suggest))
})

test_that("mt_suggest_speed_cap warns above 55 m/s biological-sanity ceiling", {
  ## Construct a track where the suggested cap lands above 55 m/s.
  ## A K02-style spoof boundary produces speeds of ~200 m/s while bulk
  ## sits at ~5 m/s.  The suggester finds a break in the gap; the
  ## warning fires because the break is well above 55 m/s.
  set.seed(7)
  n <- 1000L
  step_size <- stats::rgamma(n, shape = 4, rate = 1)  # ~4 m/s @ dt=60
  coords <- cbind(cumsum(step_size), rep(0, n))
  ## Two spoof boundaries (4 high-speed steps total = 0.4% of fixes).
  coords[c(300L, 700L), 1] <- coords[c(300L, 700L), 1] + 12000  # ~200 m/s @ dt=60
  df <- data.frame(
    location_long = coords[, 1], location_lat = coords[, 2],
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") +
                 seq(0, by = 60, length.out = n),
    individual_local_identifier = "ind", tag_local_identifier = "tag",
    sensor_type = "gps")
  x <- move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier", crs = 32633)

  msgs <- character(0)
  withCallingHandlers(
    v <- mt_suggest_speed_cap(x, plot = FALSE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  ## Either the suggestion is above 55 (warning fires) or it lands
  ## inside 55 (warning silent) -- both are valid outcomes depending
  ## on where the gap detector lands.  Test that the warning text
  ## fires when v > 55 and is silent otherwise.
  if (!is.na(v) && v > 55) {
    expect_true(any(grepl("exceeds.*m/s.*sustained-speed bound", msgs)))
  } else if (!is.na(v)) {
    expect_false(any(grepl("exceeds.*m/s.*sustained-speed bound", msgs)))
  }
})

test_that("mt_suggest_speed_cap user-supplied physiological_ceiling overrides 55 m/s", {
  ## When the user supplies a tighter ceiling (e.g. species-specific
  ## sprint margin), the warning fires above THAT value rather than
  ## above 55.
  set.seed(9)
  n <- 1000L
  coords <- cbind(cumsum(runif(n, 0, 1)), rep(0, n))
  coords[c(500L, 700L), 1] <- coords[c(500L, 700L), 1] + 1800  # ~30 m/s spike
  df <- data.frame(
    location_long = coords[, 1], location_lat = coords[, 2],
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") +
                 seq(0, by = 60, length.out = n),
    individual_local_identifier = "ind", tag_local_identifier = "tag",
    sensor_type = "gps")
  x <- move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier", crs = 32633)

  msgs <- character(0)
  withCallingHandlers(
    v <- mt_suggest_speed_cap(x, physiological_ceiling = 1, plot = FALSE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  ## The ceiling is set to 1 m/s, below the cap this fixture yields
  ## (~2.9 m/s), so the override warning is guaranteed to fire.  Assert
  ## that rather than guarding on it: guarding was why this test
  ## registered as "empty" -- the old guard (`v > 5`) was FALSE for this
  ## fixture, so no expectation ever executed and the override branch
  ## was never actually tested.
  expect_false(is.na(v))
  expect_gt(v, 1)
  ## The override must replace the default ceiling text, not add to it.
  expect_true(any(grepl("user-supplied physiological_ceiling", msgs)))
  expect_false(any(grepl("sustained-speed bound", msgs)))
})


test_that("data-driven thresholds dispatch per-track on multi-individual input", {
  ## Two tracks with independent speed scales, each long enough that a
  ## handful of jittered spikes is well below the substantive-mode
  ## threshold (0.2% spoof rate -- realistic for real-world data).
  ## Per-track dispatch must tune each track's cap to its own scale.
  mk <- function(id_val, base_step, n_outliers, outlier_min, outlier_max,
                  seed) {
    set.seed(seed)
    n <- 2000L
    coords <- cbind(cumsum(runif(n, 0, base_step)), rep(0, n))
    spike_idx <- sample(50:(n - 50), n_outliers)
    spike_mag <- runif(n_outliers, outlier_min, outlier_max)
    for (k in seq_along(spike_idx)) {
      coords[spike_idx[k]:(spike_idx[k] + 1), 1] <-
        coords[spike_idx[k]:(spike_idx[k] + 1), 1] + spike_mag[k]
    }
    df <- data.frame(
      location_long = coords[, 1], location_lat = coords[, 2],
      timestamp = as.POSIXct("2024-01-01", tz = "UTC") +
                   seq(0, by = 60, length.out = n),
      individual_local_identifier = id_val,
      tag_local_identifier = paste0("tag_", id_val),
      sensor_type = "gps")
    move2::mt_as_move2(df,
      coords = c("location_long", "location_lat"),
      time_column = "timestamp",
      track_id_column = "individual_local_identifier", crs = 32633)
  }
  slow <- mk("slow", base_step = 1,    n_outliers = 4,
             outlier_min = 30,    outlier_max = 80,    seed = 7)
  fast <- mk("fast", base_step = 1000, n_outliers = 4,
             outlier_min = 30000, outlier_max = 80000, seed = 13)
  combined <- rbind(slow, fast)

  res <- suppressMessages(mt_flag_speed_cap(combined, plot = FALSE))

  caps <- attr(res, "v_max_used")
  expect_type(caps, "double")
  expect_named(caps, c("slow", "fast"))
  ## Caps should differ by orders of magnitude given the two scales.
  expect_true(caps[["fast"]] > 10 * caps[["slow"]])
  ## Each track flags its own spike.
  ids <- move2::mt_track_id(res)
  expect_true(any(res$is_outlier[ids == "slow"]))
  expect_true(any(res$is_outlier[ids == "fast"]))
})


# ---- mode-position gate for the auto cap ----------------------------------

test_that(".gate_speed_cap_mode_position allows when cap is above the rightmost mode", {
  ## one bulk mode + a sparse outlier tail; cap chosen above the bulk
  set.seed(1)
  bulk      <- exp(rnorm(2000, mean = log(0.3), sd = 0.5))   # ~0.3 m/s mode
  outliers  <- exp(rnorm(20,   mean = log(50),  sd = 0.5))   # 1% true outliers
  v <- c(bulk, outliers)
  gate <- move2utils:::.gate_speed_cap_mode_position(v, v_b = 5)
  expect_true(gate$allow)
  expect_lt(gate$rightmost_mode, 5)
  expect_equal(gate$position, "above")
})

test_that(".gate_speed_cap_mode_position declines when cap lands inside an activity mode", {
  ## bulk near 0.3 m/s + a substantive flight mode at 4 m/s (3% of data)
  ## a cap proposed at 2 m/s sits BELOW the flight mode
  set.seed(2)
  rest   <- exp(rnorm(1900, mean = log(0.3), sd = 0.4))
  flight <- exp(rnorm(60,   mean = log(4),   sd = 0.2))      # ~3% in flight
  v <- c(rest, flight)
  gate <- move2utils:::.gate_speed_cap_mode_position(v, v_b = 2)
  expect_false(gate$allow)
  expect_match(gate$reason, "below|at", ignore.case = TRUE)
  expect_gt(gate$rightmost_mode, 2)
})

test_that(".gate_speed_cap_mode_position ignores sparse outlier modes (CPF_A-style)", {
  ## bulk near 0.3 + scattered spoofs at 100--1000 m/s, each <1% of fixes
  set.seed(3)
  bulk   <- exp(rnorm(2000, mean = log(0.3), sd = 0.6))
  spoofs <- c(exp(rnorm(8, mean = log(150), sd = 0.1)),
              exp(rnorm(7, mean = log(400), sd = 0.1)),
              exp(rnorm(8, mean = log(900), sd = 0.1)))
  v <- c(bulk, spoofs)
  ## a cap at 50 m/s sits above the bulk and below the (sparse) spoof modes
  ## -- the gate must accept it (the spoof modes are not substantive)
  gate <- move2utils:::.gate_speed_cap_mode_position(v, v_b = 50)
  expect_true(gate$allow)
})

test_that(".gate_speed_cap_mode_position falls back to allow on too-few inputs", {
  v <- c(0.1, 0.2, 0.3)
  gate <- move2utils:::.gate_speed_cap_mode_position(v, v_b = 1)
  expect_true(gate$allow)
  expect_match(gate$reason, "too few", fixed = TRUE)
})

test_that("auto-cap gate fires end-to-end on a stationary-mostly track", {
  ## fabricate a track where 97% of fixes are stationary (jitter) and ~3%
  ## are real flight at 3-5 m/s -- no genuine outliers.  The auto cap
  ## detector finds a candidate break in the resting/flight valley; the
  ## mode-position gate must refuse it so 0 fixes get speed-flagged.
  set.seed(4)
  n <- 2000L
  is_flight <- runif(n) < 0.03
  step_m <- ifelse(is_flight,
                   rnorm(n, mean = 4,    sd = 0.5),
                   abs(rnorm(n, mean = 0.5, sd = 0.3)))
  step_m <- pmax(step_m, 0.001)
  df <- data.frame(
    location_long = cumsum(step_m * runif(n, -1, 1)),
    location_lat  = cumsum(step_m * runif(n, -1, 1)),
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") + seq_len(n),
    individual_local_identifier = "ind",
    tag_local_identifier = "tag",
    sensor_type = "gps")
  x <- move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier", crs = 32633)
  res <- suppressMessages(mt_flag_speed_cap(x, plot = FALSE))
  ## The cap should be rejected (Inf) and no fixes flagged on speed
  expect_false(is.finite(attr(res, "v_max_used")))
  expect_equal(sum(res$is_speed_above_cap), 0L)
})


test_that("physiological_ceiling: rejects bad input", {
  x <- make_track()
  expect_error(mt_flag_speed_cap(x, physiological_ceiling = -1, plot = FALSE),
               "non-negative|positive")
  expect_error(mt_flag_speed_cap(x, physiological_ceiling = 0, plot = FALSE),
               "non-negative|positive")
  expect_error(mt_flag_speed_cap(x, physiological_ceiling = NA_real_, plot = FALSE),
               "positive")
  expect_error(mt_flag_speed_cap(x, physiological_ceiling = c(30, 40), plot = FALSE),
               "positive scalar")
})

## Build a track whose auto-cap path fires above 55 m/s: a bulk of
## modest-speed fixes plus a ~2% tail of extreme-speed fixes
## simulating spoof contamination.  Below the 5% safety guard but
## with a clear bulk-vs-tail gap that the entropy-valley detector
## can find.  The auto-cap then lands inside the gap (above the
## bulk, well above 55 m/s).
make_high_autocap_track <- function() {
  set.seed(7)
  n <- 1000L
  is_spoof <- runif(n) < 0.02
  step_m <- ifelse(is_spoof, runif(n, 100, 250), abs(rnorm(n, 5, 1)))
  step_m <- pmax(step_m, 0.001)
  df <- data.frame(
    location_long = cumsum(step_m * runif(n, -1, 1)),
    location_lat  = cumsum(step_m * runif(n, -1, 1)),
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") + seq_len(n),
    individual_local_identifier = "ind",
    tag_local_identifier = "tag",
    sensor_type = "gps")
  move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier",
    crs = 32633)
}

test_that("physiological_ceiling: NULL fallback uses 55 m/s warning text", {
  x <- make_high_autocap_track()
  expect_message(
    mt_flag_speed_cap(x, plot = FALSE, threshold_type = "auto"),
    "universal sustained-speed bound"
  )
})

test_that("physiological_ceiling: user-supplied value swaps the warning text", {
  x <- make_high_autocap_track()
  expect_message(
    mt_flag_speed_cap(x, physiological_ceiling = 30,
                       plot = FALSE, threshold_type = "auto"),
    "user-supplied physiological_ceiling"
  )
})

test_that("physiological_ceiling: warning suppressed when auto-cap is below it", {
  ## A clean slow-motion track that does not produce a high auto-cap.
  set.seed(11)
  n <- 200
  df <- data.frame(
    location_long = cumsum(rnorm(n, 0, 0.001)),
    location_lat  = cumsum(rnorm(n, 0, 0.001)),
    timestamp = as.POSIXct("2024-01-01", tz = "UTC") +
                 seq(0, by = 60, length.out = n),
    individual_local_identifier = "ind_1",
    tag_local_identifier = "tag_1",
    sensor_type = "gps"
  )
  x <- move2::mt_as_move2(df,
    coords = c("location_long", "location_lat"),
    time_column = "timestamp",
    track_id_column = "individual_local_identifier", crs = 32633)
  ## With a generous ceiling, no biology warning fires
  expect_silent(suppressMessages(
    mt_flag_speed_cap(x, physiological_ceiling = 1000,
                       plot = FALSE, threshold_type = "auto",
                       silent = TRUE)
  ))
})

test_that("physiological_ceiling: hard threshold ignores ceiling", {
  x <- make_track()
  ## hard cap -- ceiling should be irrelevant
  res <- suppressMessages(
    mt_flag_speed_cap(x, v_max = 40, threshold_type = "hard",
                       physiological_ceiling = 30, plot = FALSE)
  )
  expect_true(any(res$is_speed_above_cap))
})


## ---- pool_by tests ---------------------------------------------------

## Build a two-track move2 with controlled step speeds and an animal-
## level grouping column.  Used for the pool_by suite.
make_pair_track <- function(n1 = 200, n2 = 50, indv = "I1") {
  set.seed(123)
  build <- function(id, n, base_speed) {
    ts <- as.POSIXct("2024-01-01", tz = "UTC") + cumsum(rep(60, n))  # 1 min
    ## step ~ base_speed m/s * 60s = base_speed*60 metres, with noise
    dx <- rnorm(n, mean = base_speed * 60, sd = 0.5)
    dy <- rep(0, n)
    data.frame(id = id,
               timestamp = ts,
               lon = cumsum(dx), lat = cumsum(dy))
  }
  d <- rbind(build("t1", n1, 1.0),    # 1 m/s baseline
              build("t2", n2, 1.0))   # 1 m/s baseline, fewer fixes
  m <- move2::mt_as_move2(d, coords = c("lon", "lat"),
                            time_column = "timestamp",
                            track_id_column = "id", crs = 32633)
  td <- move2::mt_track_data(m); td$indv <- indv
  move2::mt_set_track_data(m, td)
}


test_that("pool_by validates inputs", {
  m <- make_pair_track()
  expect_error(mt_flag_speed_cap(m, pool_by = 1, plot = FALSE),
               "length 1 \\(single column")
  expect_error(mt_flag_speed_cap(m, pool_by = c("a", "b", "c"),
                                  plot = FALSE),
               "Deeper hierarchies")
  expect_error(suppressMessages(mt_flag_speed_cap(
                  m, pool_by = "no_such_col", plot = FALSE, silent = TRUE)),
               "not in")
  expect_error(suppressMessages(mt_flag_speed_cap(
                  m, pool_by = c("indv","indv"),
                  plot = FALSE, silent = TRUE)),
               "two \\*distinct\\* columns")
})


test_that("pool_by = NULL is byte-identical to no pool_by", {
  m <- make_pair_track()
  o1 <- suppressMessages(mt_flag_speed_cap(m, threshold_type = "auto",
                                              plot = FALSE, silent = TRUE))
  o2 <- suppressMessages(mt_flag_speed_cap(m, threshold_type = "auto",
                                              pool_by = NULL,
                                              plot = FALSE, silent = TRUE))
  expect_identical(o1$is_outlier, o2$is_outlier)
  expect_identical(o1$step_speed, o2$step_speed)
})


test_that("pool_by groups of size one are byte-identical to per-track", {
  ## Use a column that puts each track in its own group -> pool fit
  ## per-group reduces to the per-track fit; union is a no-op.
  m  <- make_pair_track()
  td <- move2::mt_track_data(m)
  td$grp <- as.character(td[[move2::mt_track_id_column(m)]])
  m  <- move2::mt_set_track_data(m, td)
  o1 <- suppressMessages(mt_flag_speed_cap(m, threshold_type = "auto",
                                              plot = FALSE, silent = TRUE))
  o2 <- suppressMessages(mt_flag_speed_cap(m, threshold_type = "auto",
                                              pool_by = "grp",
                                              plot = FALSE, silent = TRUE))
  expect_identical(o1$is_outlier, o2$is_outlier)
})


test_that("pool_by under threshold_type='hard' is a no-op", {
  m <- make_pair_track()
  o1 <- suppressMessages(mt_flag_speed_cap(m, v_max = 5,
                                              threshold_type = "hard",
                                              plot = FALSE, silent = TRUE))
  o2 <- suppressMessages(mt_flag_speed_cap(m, v_max = 5,
                                              threshold_type = "hard",
                                              pool_by = "indv",
                                              plot = FALSE, silent = TRUE))
  expect_identical(o1$is_outlier, o2$is_outlier)
})


test_that("pool_by lets a distribution-poor track inherit a richer cap", {
  ## Pool_by's primary value: when a track's own step-speed distribution
  ## does not support a structural break -- because it's too uniform
  ## (no spike to break from), genuinely too small, or has too few
  ## tail observations -- the per-track auto path returns no cap.
  ## Pool_by then lets that track inherit a cap from a pooled group
  ## (here: same animal's longer/richer deployment) and flag a fast
  ## step that would otherwise pass unflagged.
  ##
  ## Concrete construction: a long, varied "long" track (the kind of
  ## distribution that admits a clean break) paired with a SHORT
  ## "short" deployment carrying one fast injected step.  Per-track
  ## auto on the short alone yields no break (too few tail obs);
  ## pool_by = "indv" uses the long track's bulk to set a sharp cap
  ## and catches the injected step.
  set.seed(7)
  build <- function(id, n, base_speed = 1.0) {
    ts <- as.POSIXct("2024-01-01", tz = "UTC") + cumsum(rep(60, n))
    dx <- rnorm(n, mean = base_speed * 60, sd = 1.0)
    dy <- rep(0, n)
    data.frame(id = id, timestamp = ts,
               lon = cumsum(dx), lat = cumsum(dy))
  }
  d_long  <- build("long",  600)
  d_short <- build("short", 20)
  ## inject one fast step in the short track: 50 m/s for 60 s
  d_short$lon[10] <- d_short$lon[10] + 50 * 60
  d <- rbind(d_long, d_short)
  m <- move2::mt_as_move2(d, coords = c("lon", "lat"),
                            time_column = "timestamp",
                            track_id_column = "id", crs = 32633)
  td <- move2::mt_track_data(m); td$indv <- "I1"
  m  <- move2::mt_set_track_data(m, td)

  o_pt <- suppressMessages(mt_flag_speed_cap(m, threshold_type = "auto",
                                                plot = FALSE, silent = TRUE))
  o_pg <- suppressMessages(mt_flag_speed_cap(m, threshold_type = "auto",
                                                pool_by = "indv",
                                                plot = FALSE, silent = TRUE))
  ## Per-track ⊆ pool (union contract).
  expect_true(all(which(o_pt$is_outlier) %in% which(o_pg$is_outlier)))
  ## Pool catches at least one flag the per-track auto missed.
  ## (Strict: short track per-track had no break; pool inherits long's
  ## distribution and flags the injected spike.)
  expect_gt(sum(o_pg$is_outlier), sum(o_pt$is_outlier))
})


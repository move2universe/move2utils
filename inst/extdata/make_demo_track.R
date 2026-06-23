## ===========================================================================
## make_demo_track() -- a reproducible, fully synthetic GPS track for the
## move2utils outlier-cleaning vignettes.  No data file is shipped; the
## vignettes source() this script and call make_demo_track().
##
## The track is built in three layers, mirroring real telemetry:
##   1. TRUE MOVEMENT  -- a multi-state correlated random walk (rest /
##      forage / commute), so step speed is genuinely multimodal.
##   2. GPS NOISE      -- isotropic Gaussian jitter on every fix, with the
##      1-sigma scaled by fix quality (satellite count).  This is the
##      "noise floor" the detectors must leave alone.
##   3. INJECTED ERRORS -- one of each gross error class move2utils
##      distinguishes, recorded in a ground-truth column so the vignettes
##      can show exactly what was caught and what was (knowably) missed.
##
## The animal is fictional; speeds are illustrative, not a real species.
## Returns a single-individual move2 object (lon/lat) with columns:
##   gps_satellite_count, true_state, error_type ("clean" or the class).
## ===========================================================================

make_demo_track <- function(seed = 1) {
  set.seed(seed)

  ## ---- 1. irregular time schedule (15-min base, a gap, a burst) ----------
  n      <- 1500L
  t0     <- as.POSIXct("2024-05-01 04:00:00", tz = "UTC")
  dt_min <- rep(15, n - 1L)
  dt_min[600:611]  <- 1                       # a 1-min burst (12 fixes)
  dt_min[900]      <- 12 * 60                  # a ~12 h overnight gap
  dt_min[sample(seq_len(n - 1L), 60)] <- 30    # scattered missed fixes
  time <- t0 + cumsum(c(0, dt_min)) * 60
  dt_s <- c(diff(as.numeric(time)), NA)        # per-fix lag (s)

  ## ---- 2. multi-state movement (rest / forage / commute) -----------------
  ## state-specific cruising speed (m/s) and angular persistence (kappa).
  ## Speeds are modest so that normal steps stay well below the gross
  ## errors injected later (commute ~2.5 m/s -> ~2.2 km per 15-min step).
  sp_state <- c(rest = 0.02, forage = 0.5, commute = 2.5)
  ka_state <- c(rest = 0,    forage = 0.6, commute = 6)     # von-Mises-ish
  trans <- rbind(rest    = c(0.94, 0.05, 0.01),
                 forage  = c(0.03, 0.93, 0.04),
                 commute = c(0.01, 0.07, 0.92))
  state <- integer(n); state[1] <- 2L
  for (i in 2:n) state[i] <- sample.int(3L, 1L, prob = trans[state[i - 1L], ])
  state_name <- factor(c("rest", "forage", "commute")[state],
                       levels = c("rest", "forage", "commute"))

  ## State-switching correlated walk with three ingredients:
  ##  - HOME-RANGE reversion (OU pull to a centre) so the animal stays
  ##    bounded rather than diffusing away;
  ##  - GAMMA step speed per state (most steps short, occasionally long);
  ##  - VELOCITY persistence (consecutive speeds correlated) and DIRECTIONAL
  ##    persistence (correlated heading) -- so a long step sits among other
  ##    long steps and does not read as an isolated outlier.
  home_x <- 0; home_y <- 0
  beta    <- 1.1e-4                             # home-range reversion (1/s)
  rho_v   <- 0.85                               # velocity (speed) persistence
  k_shape <- 3                                  # Gamma shape (right-skew)
  heading <- numeric(n); heading[1] <- runif(1, 0, 2 * pi)
  speed   <- numeric(n); speed[1]   <- sp_state[state[1]]
  x <- numeric(n); y <- numeric(n)             # TRUE positions, metres
  for (i in 2:n) {
    ka       <- ka_state[state[i]]
    turn     <- if (ka > 0) rnorm(1, 0, 1 / sqrt(ka)) else runif(1, -pi, pi)
    heading[i] <- heading[i - 1L] + turn
    draw     <- rgamma(1, shape = k_shape, rate = k_shape / sp_state[state[i]])
    speed[i] <- rho_v * speed[i - 1L] + (1 - rho_v) * draw      # persistent m/s
    step     <- speed[i] * dt_s[i - 1L]
    x[i] <- x[i - 1L] + step * cos(heading[i]) - beta * dt_s[i - 1L] * (x[i - 1L] - home_x)
    y[i] <- y[i - 1L] + step * sin(heading[i]) - beta * dt_s[i - 1L] * (y[i - 1L] - home_y)
  }

  ## ---- 3. GPS noise, scaled by fix quality -------------------------------
  ## satellite count: mostly good, some marginal, a few poor.
  sats  <- sample(c(11, 10, 9, 8, 7, 6, 5, 4), n, replace = TRUE,
                  prob = c(.18, .2, .2, .15, .1, .07, .06, .04))
  sigma <- pmin(45, 6 + 120 / sats^1.4)        # 1-sigma metres, ~8-30 m
  xo <- x + rnorm(n, 0, sigma)                 # OBSERVED positions
  yo <- y + rnorm(n, 0, sigma)

  ## ---- 4. inject the gross-error catalogue -------------------------------
  ## Magnitudes are large relative to normal steps (commute ~2.2 km), so the
  ## gross errors are unambiguous; the subtle and halo classes sit close to
  ## the noise floor on purpose.
  err <- rep("clean", n)
  bump <- function(i, dx, dy, label) {
    xo[i] <<- xo[i] + dx; yo[i] <<- yo[i] + dy; err[i] <<- label
  }
  ## (a) out-and-back excursion: one fix darts ~9 km out and the path
  ##     returns -> path >> displacement (the detour detector's target).
  bump(330L, 8000, 4000, "out_and_back")
  ## (b) teleports / speed jumps: single impossible-speed fixes (~100 km).
  bump(455L,  90000, -60000, "teleport")
  bump(1180L, -80000, 100000, "teleport")
  ## (c) coherent SPOOF block: a run of consecutive fixes all shifted by
  ##     the SAME large offset (~46 km) -> internally consistent, only the
  ##     boundary steps look fast (the block-expansion / per-fix-evades case).
  sp <- 700:709
  xo[sp] <- xo[sp] + 30000; yo[sp] <- yo[sp] - 35000; err[sp] <- "spoof"
  ## (d) GPS-JAMMING block: a burst frozen ~50 km off the true position.
  jm <- 980:991
  xo[jm] <- x[980L] + 40000 + rnorm(length(jm), 0, 30)
  yo[jm] <- y[980L] - 30000 + rnorm(length(jm), 0, 30)
  err[jm] <- "jam"
  ## (e) COLONY HALO (perpendicular drift): during a foraging run near a
  ##     structure, multipath pushes each fix SIDEWAYS -- perpendicular to
  ##     the local travel direction -- by ~0.6-1.2 km.  Because the heading
  ##     varies through the run, the offsets point in many directions and
  ##     form a ring/halo; because each offset is orthogonal to travel, the
  ##     DIRECTIONAL bridge (dBGB, eta_perp) is the detector that sees it.
  for_runs <- rle(state == 2L)
  ends   <- cumsum(for_runs$lengths)
  starts <- ends - for_runs$lengths + 1L
  pick   <- which(for_runs$values & for_runs$lengths >= 8L)[1]
  if (!is.na(pick)) {
    hl   <- starts[pick]:(starts[pick] + 7L)
    perp <- heading[hl] + pi / 2               # orthogonal to travel
    off  <- runif(length(hl), 600, 1200) * sample(c(-1, 1), length(hl), TRUE)
    xo[hl] <- xo[hl] + off * cos(perp)
    yo[hl] <- yo[hl] + off * sin(perp)
    err[hl] <- "halo"
  }
  ## (f) subtle in-bulk errors: displaced only ~4-5 sigma -> borderline,
  ##     some genuinely unknowable from jitter (the honest-limits point).
  sub <- c(220L, 1320L)
  xo[sub] <- xo[sub] + 4.5 * sigma[sub] * sample(c(-1, 1), length(sub), TRUE)
  yo[sub] <- yo[sub] + 4.5 * sigma[sub] * sample(c(-1, 1), length(sub), TRUE)
  err[sub] <- "subtle"

  ## ---- 5. raw-download artefacts: empty geoms + a duplicate timestamp ----
  empty_idx <- c(120L, 770L, 1410L)            # failed fixes (no location)
  time[1001L] <- time[1000L]                   # duplicate timestamp

  ## ---- 6. metres -> lon/lat and assemble a move2 -------------------------
  ref_lon <- 11.50; ref_lat <- 47.50
  lon <- ref_lon + xo / (111000 * cos(ref_lat * pi / 180))
  lat <- ref_lat + yo / 111000

  sf_obj <- sf::st_as_sf(
    data.frame(lon = lon, lat = lat,
               timestamp           = time,
               individual          = "sim_animal_01",
               gps_satellite_count = sats,
               true_state          = state_name,
               error_type          = err),
    coords = c("lon", "lat"), crs = 4326)
  ## blank the failed fixes to empty geometry
  sf::st_geometry(sf_obj)[empty_idx] <-
    sf::st_geometry(sf::st_sfc(sf::st_point(), crs = 4326))
  sf_obj$error_type[empty_idx] <- "empty_geometry"

  move2::mt_as_move2(sf_obj, time_column = "timestamp",
                     track_id_column = "individual")
}

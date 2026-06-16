## ==========================================================================
## Generate synthetic GPS tracking data for move2utils worked examples
##
## Produces Movebank-format CSV files readable by move2::mt_read():
##   synthetic_tracks.csv.gz — one file with three individuals:
##     CPF_A: central-place forager, ~2000 locations, 30 days
##            with faint, moderate, and strong isolated outliers
##            plus an outlier burst (spoofing episode)
##     CPF_B: same process, ~4000 locations, 60 days, clean
##            (serves as a reference track)
##     CPF_C: same process, ~200 locations, 3 days, short track
##            with a few moderate outliers
##
## All tracks have IRREGULAR sampling modelled after real GPS data:
## base rate ~15 min, with missed fixes, overnight gaps, and
## occasional burst sampling.
##
## Movement is simulated from a ctmm OUF model with known parameters,
## enabling comparison of recovered distributions against theory.
## ==========================================================================

library(ctmm)
set.seed(2025)

## ---- OUF model (shared "species") ------------------------------------------
## tau_position = 2 days  (home-range crossing time)
## tau_velocity = 2 hours (directional persistence)
## sigma = 50 km^2       (position variance → HR area ≈ 941 km^2)
## Characteristic speed ≈ 24.5 km/day
ouf_model <- ctmm(tau = c(48 %#% "hour", 2 %#% "hour"),
                   sigma = 50 %#% "km^2",
                   mu = c(0, 0),
                   isotropic = TRUE)

## geographic reference (Alps region)
ref_lon <- 11.50
ref_lat <- 47.50
m_per_deg_lon <- 111000 * cos(ref_lat * pi / 180)
m_per_deg_lat <- 111000

cat("OUF model:\n")
print(summary(ouf_model))


## ---- helper: generate irregular time schedule ------------------------------
## Mimics real GPS data: base rate with missed fixes, overnight gaps,
## and occasional burst sampling
make_schedule <- function(n_days, base_dt_min = 15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  ## generate regular schedule, then thin and add jitter
  n_regular <- n_days * 24 * 60 / base_dt_min
  t_regular <- seq(0, by = base_dt_min * 60, length.out = n_regular)

  ## keep each fix with probability depending on time of day
  ## (lower at night = duty cycling / no solar power)
  hours <- (t_regular / 3600) %% 24
  p_keep <- ifelse(hours >= 6 & hours <= 20,
                   0.85,   # daytime: 85% fix rate
                   0.30)   # nighttime: 30% fix rate

  ## add some random dropout clusters (simulates poor satellite geometry)
  n_dropout <- rpois(1, n_days / 3)
  if (n_dropout > 0) {
    dropout_centres <- sample(length(t_regular), n_dropout)
    for (dc in dropout_centres) {
      window <- max(1, dc - 8):min(length(t_regular), dc + 8)
      p_keep[window] <- p_keep[window] * 0.2
    }
  }

  keep <- runif(length(t_regular)) < p_keep
  keep[1] <- TRUE  # always keep the first fix
  t_sched <- t_regular[keep]

  ## add small timestamp jitter (±30 seconds, realistic for GPS)
  t_sched <- t_sched + runif(length(t_sched), -30, 30)
  t_sched <- sort(t_sched)
  t_sched[1] <- 0  # anchor start

  t_sched
}


## ---- helper: convert ctmm simulation to lon/lat ----------------------------
ctmm_to_lonlat <- function(sim) {
  data.frame(
    lon = ref_lon + sim$x / m_per_deg_lon,
    lat = ref_lat + sim$y / m_per_deg_lat
  )
}


## ---- helper: build Movebank-format data frame ------------------------------
make_movebank_df <- function(coords, timestamps, id,
                             sat_count = NULL, dop = NULL) {
  n <- nrow(coords)
  if (is.null(sat_count)) sat_count <- sample(5:12, n, replace = TRUE)
  if (is.null(dop))       dop <- round(runif(n, 1.0, 4.0), 1)

  ## ground speed (m/s) from consecutive positions (approximate)
  dx_m <- c(NA, diff(coords$lon)) * m_per_deg_lon
  dy_m <- c(NA, diff(coords$lat)) * m_per_deg_lat
  dist_m <- sqrt(dx_m^2 + dy_m^2)
  dt_s <- c(NA, as.numeric(diff(timestamps), units = "secs"))
  gspeed <- ifelse(dt_s > 0, dist_m / dt_s, 0)
  gspeed[1] <- 0

  hdg <- atan2(dx_m, dy_m) * 180 / pi
  hdg <- ifelse(hdg < 0, hdg + 360, hdg)
  hdg[1] <- 0

  data.frame(
    timestamp = format(timestamps, "%Y-%m-%d %H:%M:%S.000"),
    `location-long` = coords$lon,
    `location-lat`  = coords$lat,
    `ground-speed`  = round(gspeed, 2),
    heading       = round(hdg, 0),
    `gps:satellite-count` = sat_count,
    `gps:dop`       = dop,
    `individual-local-identifier` = id,
    `sensor-type`   = "gps",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


## ---- helper: inject outliers -----------------------------------------------
## Returns modified coords and a ground-truth data frame
inject_outliers <- function(coords, spec) {
  ## spec is a list of outlier specifications, each with:
  ##   indices: which positions to displace
  ##   magnitude_deg: displacement in degrees (approximate)
  ##   type: "faint", "moderate", "strong", "burst"
  gt <- data.frame(index = integer(), type = character(),
                   magnitude = numeric(), stringsAsFactors = FALSE)

  for (s in spec) {
    for (idx in s$indices) {
      angle <- runif(1, 0, 2 * pi)
      mag <- s$magnitude_deg * runif(1, 0.7, 1.3)  # ±30% variation
      coords$lon[idx] <- coords$lon[idx] + mag * cos(angle)
      coords$lat[idx] <- coords$lat[idx] + mag * sin(angle)
      gt <- rbind(gt, data.frame(index = idx, type = s$type,
                                  magnitude = mag))
    }
  }
  list(coords = coords, ground_truth = gt)
}


## ============================================================================
## Track A: central-place forager with diverse outlier types
## ============================================================================
cat("\n--- Track A: CPF with diverse outliers ---\n")

t_start_a <- as.POSIXct("2025-06-01 06:00:00", tz = "UTC")
sched_a <- make_schedule(30, base_dt_min = 15, seed = 101)
t_seq_a <- as.numeric(t_start_a) + sched_a

cat("  Schedule: ", length(t_seq_a), " locations over 30 days\n")
cat("  Sampling interval quantiles (min):\n")
print(round(quantile(diff(sched_a) / 60,
                      probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)), 1))

sim_a <- simulate(ouf_model, t = t_seq_a, seed = 42)
coords_a <- ctmm_to_lonlat(sim_a)
timestamps_a <- as.POSIXct(t_seq_a, origin = "1970-01-01", tz = "UTC")

## Normal step size for reference (median step in degrees)
step_deg <- median(sqrt(diff(coords_a$lon)^2 + diff(coords_a$lat)^2),
                   na.rm = TRUE)
cat("  Median step size:", round(step_deg, 5), "degrees\n")
cat("  ≈", round(step_deg * 111000, 0), "metres\n")

## Define outlier levels relative to the median step:
##   faint:    ~10× median step (subtle, near the edge of natural variation)
##   moderate: ~50× median step (clearly anomalous but within a few km)
##   strong:   ~500× median step (tens of km, obviously wrong)
##   burst:    ~2000× median step (hundreds of km, spoofing-like)
n_a <- length(t_seq_a)
outlier_spec_a <- list(
  ## faint isolated outliers (3 of them)
  list(indices = c(150, 600, 1200),
       magnitude_deg = step_deg * 10,
       type = "faint"),
  ## moderate isolated outliers (3 of them)
  list(indices = c(300, 900, 1500),
       magnitude_deg = step_deg * 50,
       type = "moderate"),
  ## strong isolated outliers (2 of them)
  list(indices = c(500, 1100),
       magnitude_deg = step_deg * 500,
       type = "strong"),
  ## burst: 15 consecutive locations displaced far (spoofing episode)
  list(indices = 1000:1014,
       magnitude_deg = step_deg * 2000,
       type = "burst")
)

result_a <- inject_outliers(coords_a, outlier_spec_a)
coords_a <- result_a$coords
gt_a <- result_a$ground_truth

## GPS metadata: burst locations get plausible-looking quality
sat_a <- sample(5:12, n_a, replace = TRUE)
dop_a <- round(runif(n_a, 1.0, 4.0), 1)
sat_a[1000:1014] <- sample(4:8, 15, replace = TRUE)
dop_a[1000:1014] <- round(runif(15, 2.0, 6.0), 1)

df_a <- make_movebank_df(coords_a, timestamps_a, id = "CPF_A",
                          sat_count = sat_a, dop = dop_a)

cat("  Outliers injected:\n")
cat("    faint:", sum(gt_a$type == "faint"), "\n")
cat("    moderate:", sum(gt_a$type == "moderate"), "\n")
cat("    strong:", sum(gt_a$type == "strong"), "\n")
cat("    burst:", sum(gt_a$type == "burst"), "\n")
cat("    total:", nrow(gt_a), "of", n_a, "locations\n")


## ============================================================================
## Track B: clean reference track (longer, same process)
## ============================================================================
cat("\n--- Track B: clean reference ---\n")

t_start_b <- as.POSIXct("2025-04-01 06:00:00", tz = "UTC")
sched_b <- make_schedule(60, base_dt_min = 15, seed = 202)
t_seq_b <- as.numeric(t_start_b) + sched_b

sim_b <- simulate(ouf_model, t = t_seq_b, seed = 123)
coords_b <- ctmm_to_lonlat(sim_b)
timestamps_b <- as.POSIXct(t_seq_b, origin = "1970-01-01", tz = "UTC")
df_b <- make_movebank_df(coords_b, timestamps_b, id = "CPF_B")
cat("  ", length(t_seq_b), "locations, clean\n")


## ============================================================================
## Track C: short track with moderate outliers
## ============================================================================
cat("\n--- Track C: short track with outliers ---\n")

t_start_c <- as.POSIXct("2025-07-01 06:00:00", tz = "UTC")
sched_c <- make_schedule(3, base_dt_min = 15, seed = 303)
t_seq_c <- as.numeric(t_start_c) + sched_c

sim_c <- simulate(ouf_model, t = t_seq_c, seed = 456)
coords_c <- ctmm_to_lonlat(sim_c)
timestamps_c <- as.POSIXct(t_seq_c, origin = "1970-01-01", tz = "UTC")
n_c <- length(t_seq_c)

outlier_spec_c <- list(
  list(indices = c(30, 80, 120, 160),
       magnitude_deg = step_deg * 50,
       type = "moderate")
)
result_c <- inject_outliers(coords_c, outlier_spec_c)
coords_c <- result_c$coords
gt_c <- result_c$ground_truth

df_c <- make_movebank_df(coords_c, timestamps_c, id = "CPF_C")
cat("  ", n_c, "locations, ", nrow(gt_c), "outliers\n")


## ============================================================================
## Track D: sustained spoof block (long-distance migrator + spoof event)
##
## Mimics a GNSS-spoofing event in which the receiver locks onto a false
## constellation for several hours and reports geographically displaced
## but internally coherent positions.  The block is a separate OUF
## realisation (different mu, independent stochastic state) substituted
## for an interior segment.  Expected to be caught by:
##   - bridge primitive at the spoof boundaries
##   - block expansion, which uses topological connectivity to isolate
##     the spoofed train (interior fixes look kinematically normal but
##     are disconnected from the main track through impossible jumps)
## ============================================================================
cat("\n--- Track D: spoof block on long-distance migrator ---\n")

t_start_d <- as.POSIXct("2025-09-01 06:00:00", tz = "UTC")
sched_d <- make_schedule(30, base_dt_min = 15, seed = 404)
t_seq_d <- as.numeric(t_start_d) + sched_d

## A "migrator" with longer correlation in velocity (more directed
## travel), lower position variance.
ouf_migrator <- ctmm(tau = c(120 %#% "hour", 6 %#% "hour"),
                      sigma = 100 %#% "km^2",
                      mu = c(0, 0),
                      isotropic = TRUE)
sim_d_main <- simulate(ouf_migrator, t = t_seq_d, seed = 808)
coords_d <- ctmm_to_lonlat(sim_d_main)

## Pick an interior segment: 30 consecutive fixes (~7.5 h).  Replace
## them with positions from a DIFFERENT OUF realisation centred 800 km
## NE of the main track (the "spoofed" location).
spoof_idx <- 850:879
n_spoof <- length(spoof_idx)

## Generate the spoof segment with its own model and offset
ouf_spoof <- ctmm(tau = c(48 %#% "hour", 2 %#% "hour"),
                   sigma = 5 %#% "km^2",
                   mu = c(0, 0),
                   isotropic = TRUE)
sim_d_spoof <- simulate(ouf_spoof, t = t_seq_d[spoof_idx], seed = 909)

## Offset the spoof segment 800 km NE of the main track at the time of
## the first spoofed fix
spoof_offset_x <- sim_d_main$x[spoof_idx[1]] + 800e3 * cos(pi / 4)
spoof_offset_y <- sim_d_main$y[spoof_idx[1]] + 800e3 * sin(pi / 4)
coords_d$lon[spoof_idx] <- ref_lon + (sim_d_spoof$x + spoof_offset_x) / m_per_deg_lon
coords_d$lat[spoof_idx] <- ref_lat + (sim_d_spoof$y + spoof_offset_y) / m_per_deg_lat

timestamps_d <- as.POSIXct(t_seq_d, origin = "1970-01-01", tz = "UTC")
df_d <- make_movebank_df(coords_d, timestamps_d, id = "CPF_D")
gt_d <- data.frame(index = spoof_idx, type = "spoof_block",
                    magnitude = NA_real_, stringsAsFactors = FALSE)
cat("  ", length(t_seq_d), "locations,", n_spoof, "spoof-block fixes (",
    spoof_idx[1], "-", spoof_idx[n_spoof], ")\n", sep = "")

## ============================================================================
## Track E: colony with stationary GPS-jitter halo
##
## Mimics a central-place forager at a colony with sparse 1-h GPS
## sampling and a population of injected radial-spike fixes.  At 1-h
## sampling, a 50-100 km radial spike implies step speed of 14-28 m/s
## -- below the gull physiological cap (~36 m/s) -- so the speed-cap
## detector is structurally insensitive.  Bridge perpendicular residual
## and detour ratio are both designed to catch this class.
##
## The injected spikes are placed at random throughout the track
## (typical real-world stationary GPS jitter) with magnitudes drawn
## uniformly from 50-100 km radius.
## ============================================================================
cat("\n--- Track E: colony with GPS-jitter halo (1-h sampling) ---\n")

t_start_e <- as.POSIXct("2025-05-15 06:00:00", tz = "UTC")
n_days_e <- 60
## 1-h sampling, no irregular dropout (assume colony has good signal)
sched_e <- seq(0, by = 3600, length.out = n_days_e * 24)
t_seq_e <- as.numeric(t_start_e) + sched_e

## Stationary OUF: very tight position variance ~ colony residence
ouf_stationary <- ctmm(tau = c(48 %#% "hour", 0.5 %#% "hour"),
                        sigma = 0.1 %#% "km^2",      # tight HR ~ 1.9 km^2
                        mu = c(0, 0),
                        isotropic = TRUE)
sim_e <- simulate(ouf_stationary, t = t_seq_e, seed = 505)
coords_e <- ctmm_to_lonlat(sim_e)

## Inject 80 isolated spike fixes uniformly throughout the track at
## radii 50-100 km.  Avoid the first/last 5 fixes so the injection
## doesn't sit at a boundary where bridge can't compute a residual.
n_e <- length(t_seq_e)
n_spikes <- 80
spike_idx <- sort(sample(6:(n_e - 5), n_spikes, replace = FALSE))
spike_radii_km <- runif(n_spikes, 50, 100)
spike_angles  <- runif(n_spikes, 0, 2 * pi)
for (k in seq_along(spike_idx)) {
  i <- spike_idx[k]
  r <- spike_radii_km[k] * 1000   # m
  coords_e$lon[i] <- coords_e$lon[i] + r * cos(spike_angles[k]) / m_per_deg_lon
  coords_e$lat[i] <- coords_e$lat[i] + r * sin(spike_angles[k]) / m_per_deg_lat
}
timestamps_e <- as.POSIXct(t_seq_e, origin = "1970-01-01", tz = "UTC")
df_e <- make_movebank_df(coords_e, timestamps_e, id = "CPF_E")
gt_e <- data.frame(index = spike_idx, type = "halo_spike",
                    magnitude = spike_radii_km * 1000,
                    stringsAsFactors = FALSE)
cat("  ", n_e, "locations,", n_spikes, "halo spikes (50-100 km radii)\n",
    sep = "")

## ============================================================================
## Track F: multi-state migrator (rest + flight) with state-relative
## anomalies
##
## Two segments concatenated: a stationary period (rest at colony,
## small sigma, low tau_velocity) followed by a directed migration
## (long tau_velocity, large sigma).  Within each, two outliers are
## injected:
##   - rest period: a 30-km step (anomalous for rest, normal for flight)
##   - flight period: a near-zero step (anomalous for flight, normal
##     for rest)
##
## Single-state cleaning would see a bimodal kinematic distribution and
## either over-flag the rest mode or under-flag the flight mode;
## state-conditional cleaning, given a per-fix state vector, correctly
## flags both within their own state context.
## ============================================================================
cat("\n--- Track F: multi-state migrator ---\n")

t_start_f <- as.POSIXct("2025-04-01 06:00:00", tz = "UTC")
sched_f_rest    <- seq(0, by = 900, length.out = 15 * 24 * 4)   # 15 d
sched_f_flight  <- seq(max(sched_f_rest) + 900, by = 900,
                        length.out = 15 * 24 * 4)               # 15 d
sched_f <- c(sched_f_rest, sched_f_flight)
t_seq_f <- as.numeric(t_start_f) + sched_f

## Rest segment: very stationary
sim_f_rest <- simulate(ouf_stationary, t = t_seq_f[seq_along(sched_f_rest)],
                        seed = 606)
## Flight segment: directed migration (long velocity tau, large sigma)
ouf_flight <- ctmm(tau = c(240 %#% "hour", 12 %#% "hour"),
                    sigma = 500 %#% "km^2",
                    mu = c(0, 0),
                    isotropic = TRUE)
sim_f_flight <- simulate(ouf_flight,
                          t = t_seq_f[length(sched_f_rest) +
                                      seq_along(sched_f_flight)],
                          seed = 707)
## Offset the flight segment to start where the rest segment ends so
## the concatenation is continuous in space
flight_offset_x <- sim_f_rest$x[length(sched_f_rest)] -
                   sim_f_flight$x[1]
flight_offset_y <- sim_f_rest$y[length(sched_f_rest)] -
                   sim_f_flight$y[1]
sim_f <- list(
  x = c(sim_f_rest$x, sim_f_flight$x + flight_offset_x),
  y = c(sim_f_rest$y, sim_f_flight$y + flight_offset_y)
)
coords_f <- ctmm_to_lonlat(sim_f)
n_f <- length(t_seq_f)

## State labels (per-fix)
state_f <- c(rep("rest",   length(sched_f_rest)),
             rep("flight", length(sched_f_flight)))

## State-relative outliers:
##   - rest: indices 200, 600 -- a 2 km displacement
##           (~10x rest's ~100 m typical step, anomalous in-state;
##            within flight's ~5-10 km typical, normal globally;
##            implied speed at 15-min sampling = 2.2 m/s, well
##            below the 36 m/s physiological cap so pre-peel
##            cannot catch it)
##   - flight: indices 2000, 2500 -- a sharp 90-degree turn
##           combined with reduced step magnitude (anomalous for
##           the directed-flight kinematics; normal for rest).
rest_outlier_idx   <- c(200, 600)
flight_outlier_idx <- c(2000, 2500)
for (i in rest_outlier_idx) {
  coords_f$lon[i] <- coords_f$lon[i] + 2000 / m_per_deg_lon
  coords_f$lat[i] <- coords_f$lat[i] + 2000 / m_per_deg_lat
}
## For flight outliers, displace orthogonally to the local direction by
## ~1 km, well below typical flight step but anomalous given the
## directional persistence the flight state implies.
for (i in flight_outlier_idx) {
  ## local direction at fix i-1 -> i
  dx <- coords_f$lon[i] - coords_f$lon[i - 1]
  dy <- coords_f$lat[i] - coords_f$lat[i - 1]
  ## perpendicular unit vector
  nrm <- sqrt(dx^2 + dy^2)
  if (is.finite(nrm) && nrm > 0) {
    perp_x <- -dy / nrm
    perp_y <-  dx / nrm
    coords_f$lon[i] <- coords_f$lon[i] + 1000 * perp_x / m_per_deg_lon
    coords_f$lat[i] <- coords_f$lat[i] + 1000 * perp_y / m_per_deg_lat
  }
}
timestamps_f <- as.POSIXct(t_seq_f, origin = "1970-01-01", tz = "UTC")
df_f <- make_movebank_df(coords_f, timestamps_f, id = "CPF_F")
df_f$state <- state_f       # state column travels with the data
gt_f <- data.frame(
  index     = c(rest_outlier_idx, flight_outlier_idx),
  type      = c(rep("rest_state_anomaly",   length(rest_outlier_idx)),
                rep("flight_state_anomaly", length(flight_outlier_idx))),
  magnitude = c(rep(2000, length(rest_outlier_idx)),
                rep(1000, length(flight_outlier_idx))),
  stringsAsFactors = FALSE)
cat("  ", n_f, "locations,",
    nrow(gt_f), "state-anomalous fixes (",
    length(rest_outlier_idx), " rest +",
    length(flight_outlier_idx), " flight)\n", sep = "")

## ============================================================================
## Write output
## ============================================================================
cat("\nWriting files...\n")

## combine all tracks into one file (Movebank convention)
## Tracks D/E/F may have a `state` column; A/B/C don't.  Reconcile by
## adding the column where missing.
for (df in list("df_a","df_b","df_c","df_d","df_e")) {
  d <- get(df)
  if (is.null(d$state)) d$state <- NA_character_
  assign(df, d)
}
df_all <- rbind(df_a, df_b, df_c, df_d, df_e, df_f)

outpath <- "inst/extdata/synthetic_tracks.csv.gz"
con <- gzfile(outpath, "w")
writeLines(paste(names(df_all), collapse = ","), con)
write.table(df_all, con, sep = ",", row.names = FALSE,
            col.names = FALSE, quote = TRUE)
close(con)
cat("  Written:", outpath, "\n")

## ground truth (per-track outlier indices + types)
gt <- list(
  CPF_A = gt_a,
  CPF_C = gt_c,
  CPF_D = gt_d,
  CPF_E = gt_e,
  CPF_F = gt_f,
  ouf_model = ouf_model,
  ref_point = c(lon = ref_lon, lat = ref_lat)
)
saveRDS(gt, "inst/extdata/synthetic_ground_truth.rds")
cat("  Written: synthetic_ground_truth.rds\n")

cat("\nDone. Total:", nrow(df_all), "locations across 6 tracks.\n")
cat("  CPF_A:", nrow(df_a), "(", nrow(gt_a), "outliers, mixed types)\n")
cat("  CPF_B:", nrow(df_b), "(clean reference)\n")
cat("  CPF_C:", nrow(df_c), "(", nrow(gt_c), "outliers, moderate)\n")
cat("  CPF_D:", nrow(df_d), "(", nrow(gt_d), "fixes in spoof block)\n")
cat("  CPF_E:", nrow(df_e), "(", nrow(gt_e), "GPS-jitter halo spikes)\n")
cat("  CPF_F:", nrow(df_f), "(", nrow(gt_f), "state-anomalous outliers)\n")

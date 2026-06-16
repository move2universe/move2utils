# ---- helper: simulate a correlated random walk as move2 ----

#' @noRd
.simulate_crw_move2 <- function(n = 500, step_sd = 0.001, turn_sd = 0.5,
                                 start = c(10, 48), seed = 42) {
  set.seed(seed)
  ## correlated turning angles (centred at 0 = forward)
  turns <- cumsum(rnorm(n, mean = 0, sd = turn_sd))
  ## step_sd is in degrees (approx 0.001 deg ~ 100m)
  steps <- abs(rnorm(n, mean = step_sd, sd = step_sd * 0.3))

  lon <- numeric(n + 1)
  lat <- numeric(n + 1)
  lon[1] <- start[1]
  lat[1] <- start[2]

  for (i in seq_len(n)) {
    lon[i + 1] <- lon[i] + steps[i] * cos(turns[i])
    lat[i + 1] <- lat[i] + steps[i] * sin(turns[i])
  }

  ## build move2 object (WGS84 for mt_turnangle compatibility)
  timestamps <- as.POSIXct("2024-01-01", tz = "UTC") + seq(0, n) * 3600
  df <- data.frame(lon = lon, lat = lat, timestamp = timestamps,
                   track_id = "sim1")
  pts <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  move2::mt_as_move2(pts, time_column = "timestamp",
                     track_id_column = "track_id")
}

#' Inject known outliers by displacing locations
#' @noRd
.inject_outliers <- function(mv2, n_outliers = 10, displacement = 0.05,
                              seed = 123) {
  set.seed(seed)
  n <- nrow(mv2)
  ## avoid first and last 2 locations (no turn angle there)
  candidates <- 3:(n - 1)
  inj_idx <- sort(sample(candidates, min(n_outliers, length(candidates))))

  geom <- sf::st_geometry(mv2)
  for (i in inj_idx) {
    pt <- sf::st_coordinates(geom[i])
    angle <- runif(1, 0, 2 * pi)
    new_pt <- c(pt[1] + displacement * cos(angle),
                pt[2] + displacement * sin(angle))
    geom[i] <- sf::st_point(new_pt)
  }
  sf::st_geometry(mv2) <- geom
  list(data = mv2, injected_idx = inj_idx)
}

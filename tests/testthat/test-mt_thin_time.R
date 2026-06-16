test_that("mt_thin_time rejects non-move2 input", {
  expect_error(mt_thin_time(data.frame(t = 1:3), interval = 1),
               "inherits")
})

test_that("mt_thin_time rejects non-positive interval or invalid tolerance", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M1", ][1:30, ]

  expect_error(mt_thin_time(x, interval = 0), "must be positive")
  expect_error(mt_thin_time(x, interval = 60, tolerance = -1),
               "non-negative")
  expect_error(mt_thin_time(x, interval = 60, tolerance = 120),
               "must not exceed")
})

test_that("mt_thin_time accepts a units object and converts it (not strips it)", {
  skip_if_not_installed("move2")
  skip_if_not_installed("units")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M1", ][1:100, ]

  ## a units object in minutes must thin identically to the equivalent
  ## string/difftime; the old bug stripped "45 min" down to 45 seconds.
  by_units  <- mt_thin_time(x, interval = units::set_units(45, "min"),
                            tolerance = units::set_units(5, "min"))
  by_string <- mt_thin_time(x, interval = "45 min", tolerance = "5 min")
  by_dt     <- mt_thin_time(x, interval = as.difftime(45, units = "mins"),
                            tolerance = as.difftime(5, units = "mins"))
  expect_equal(by_units$thin_selected, by_string$thin_selected)
  expect_equal(by_units$thin_selected, by_dt$thin_selected)

  ## a non-time unit cannot be coerced to seconds -> classed error
  expect_error(
    mt_thin_time(x, interval = units::set_units(45, "m")),
    class = "move2utils_mt_thin_time_bad_units")
})

test_that("mt_thin_time returns a thin_selected column of correct length", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M1", ][1:100, ]

  out <- mt_thin_time(x, interval = "10 min", tolerance = "2 min")
  expect_s3_class(out, "move2")
  expect_true("thin_selected" %in% names(out))
  expect_type(out$thin_selected, "logical")
  expect_equal(length(out$thin_selected), nrow(x))
  expect_equal(nrow(out), nrow(x))
})

test_that("mt_thin_time with remove = TRUE returns only selected rows", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M1", ][1:100, ]

  full    <- mt_thin_time(x, interval = "10 min", tolerance = "2 min")
  trimmed <- mt_thin_time(x, interval = "10 min", tolerance = "2 min",
                          remove = TRUE)
  expect_equal(nrow(trimmed), sum(full$thin_selected))
})

test_that("retained fixes never have lags shorter than interval - tolerance", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M4", ][1:200, ]

  out <- mt_thin_time(x, interval = "5 min", tolerance = "1 min",
                      remove = TRUE)
  lags <- as.numeric(move2::mt_time_lags(out, units = "min"))
  lags <- lags[!is.na(lags)]
  ## lags below the lower tolerance bound are forbidden by construction;
  ## lags above the upper bound occur only at run boundaries
  expect_true(all(lags >= 4))
})

test_that("run-splitting prevents bridging across large gaps", {
  skip_if_not_installed("move2")
  ## synthetic: two sub-runs separated by a large gap
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  times <- c(t0 + seq(0, 600, by = 60),
             t0 + 7200 + seq(0, 600, by = 60))   # 2h gap
  coords <- data.frame(x = seq_along(times),
                        y = seq_along(times))
  sf_pts <- sf::st_as_sf(coords, coords = c("x", "y"), crs = 4326)
  sf_pts$time <- times
  sf_pts$tid  <- "synthetic"
  x <- move2::mt_as_move2(sf_pts, time_column = "time",
                           track_id_column = "tid")

  out <- mt_thin_time(x, interval = "60 secs", tolerance = "10 secs",
                      remove = TRUE)
  lags <- as.numeric(move2::mt_time_lags(out, units = "secs"))
  lags <- lags[!is.na(lags)]
  ## no retained lag should land in the forbidden interval
  ## (> 70 secs and < 6600 secs — i.e. bridging across the 2h gap)
  expect_false(any(lags > 70 & lags < 6600))
})

test_that("multi-track input processes each track independently", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  ## stratified subset — take the first 100 of each individual
  subset_id <- function(id) fishers[move2::mt_track_id(fishers) == id, ][1:100, ]
  two <- rbind(subset_id("F1"), subset_id("F2"))

  out <- mt_thin_time(two, interval = "10 min", tolerance = "2 min")
  expect_true("thin_selected" %in% names(out))

  by_id <- tapply(out$thin_selected,
                  droplevels(factor(move2::mt_track_id(out))),
                  sum)
  ## each present individual should have at least one selected fix
  expect_true(all(by_id > 0))
})

test_that("criterion = 'first' vs 'closest' can differ", {
  skip_if_not_installed("move2")
  ## synthetic track where two chains of equal length exist
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  times <- t0 + c(0, 40, 50, 90, 95, 135) * 60
  coords <- data.frame(x = seq_along(times), y = seq_along(times))
  sf_pts <- sf::st_as_sf(coords, coords = c("x", "y"), crs = 4326)
  sf_pts$time <- times
  sf_pts$tid  <- "tie"
  x <- move2::mt_as_move2(sf_pts, time_column = "time",
                           track_id_column = "tid")

  a <- mt_thin_time(x, interval = "45 min", tolerance = "5 min",
                    criterion = "closest", remove = TRUE)
  b <- mt_thin_time(x, interval = "45 min", tolerance = "5 min",
                    criterion = "first",   remove = TRUE)

  ## both return valid chains; "closest" must be at least as good
  la <- abs(as.numeric(move2::mt_time_lags(a, units = "min")) - 45)
  lb <- abs(as.numeric(move2::mt_time_lags(b, units = "min")) - 45)
  expect_true(sum(la, na.rm = TRUE) <= sum(lb, na.rm = TRUE) + 1e-9)
})

## tests for mt_filter_gps_quality

make_track <- function(sat = NULL, dop = NULL, hacc = NULL, n = 20) {
  coords <- cbind(runif(n, 10, 11), runif(n, 47, 48))
  df <- data.frame(
    location_long = coords[, 1],
    location_lat  = coords[, 2],
    timestamp = as.POSIXct("2024-01-01 00:00:00", tz = "UTC") +
                 seq(0, by = 600, length.out = n),
    individual_local_identifier = "ind_1",
    tag_local_identifier = "tag_1",
    sensor_type = "gps"
  )
  if (!is.null(sat))  df$gps_satellite_count <- sat
  if (!is.null(dop))  df$gps_dop <- dop
  if (!is.null(hacc)) df$eobs_horizontal_accuracy_estimate <- hacc
  move2::mt_as_move2(df,
                     coords = c("location_long", "location_lat"),
                     time_column = "timestamp",
                     track_id_column = "individual_local_identifier",
                     crs = 4326)
}

test_that("rejects non-move2 input", {
  expect_error(mt_filter_gps_quality(data.frame(x = 1)),
               class = "move2utils_input_not_move2")
})

test_that("no-op when no quality columns present", {
  x <- make_track(n = 10)
  expect_message(y <- mt_filter_gps_quality(x), "No GPS-quality columns")
  expect_equal(nrow(y), nrow(x))
})

test_that("drops fixes below satellite threshold", {
  set.seed(1)
  sat <- c(rep(3, 5), rep(6, 15))  # five bad, fifteen good
  x <- make_track(sat = sat)
  y <- suppressMessages(mt_filter_gps_quality(x, sat_min = 5,
                                               dop_max = NULL,
                                               hacc_max = NULL))
  expect_equal(nrow(y), 15)
})

test_that("drops fixes above DOP threshold", {
  dop <- c(rep(12, 5), rep(4, 15))  # five bad
  x <- make_track(dop = dop)
  y <- suppressMessages(mt_filter_gps_quality(x, sat_min = NULL,
                                               dop_max = 10,
                                               hacc_max = NULL))
  expect_equal(nrow(y), 15)
})

test_that("drops fixes above hacc threshold", {
  hacc <- c(rep(500, 3), rep(20, 17))  # three bad
  x <- make_track(hacc = hacc)
  y <- suppressMessages(mt_filter_gps_quality(x, sat_min = NULL,
                                               dop_max = NULL,
                                               hacc_max = 100))
  expect_equal(nrow(y), 17)
})

test_that("combines criteria via AND (all must pass)", {
  ## row 1: sat=3 bad;  row 2: dop=15 bad;  row 3: hacc=500 bad;  rest good
  sat  <- c(3, 6, 6, rep(6, 17))
  dop  <- c(4, 15, 4, rep(4, 17))
  hacc <- c(20, 20, 500, rep(20, 17))
  x <- make_track(sat = sat, dop = dop, hacc = hacc)
  y <- suppressMessages(mt_filter_gps_quality(x))
  expect_equal(nrow(y), 17)
})

test_that("NULL threshold disables that criterion", {
  sat  <- c(3, rep(6, 19))   # 1 would be dropped by sat_min=5
  x <- make_track(sat = sat)
  y <- suppressMessages(mt_filter_gps_quality(x, sat_min = NULL,
                                               dop_max = NULL,
                                               hacc_max = NULL))
  expect_equal(nrow(y), 20)
})

test_that("NA values do not count as failures", {
  sat <- c(NA_real_, NA_real_, rep(6, 18))
  x <- make_track(sat = sat)
  y <- suppressMessages(mt_filter_gps_quality(x, sat_min = 5,
                                               dop_max = NULL,
                                               hacc_max = NULL))
  expect_equal(nrow(y), 20)
})

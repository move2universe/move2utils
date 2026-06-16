test_that(".bm_variance returns expected structure", {
  # Simple synthetic track: 5 locations in a line
  x <- c(0, 100, 200, 300, 400)
  y <- c(0, 0, 0, 0, 0)
  time_lag <- c(10, 10, 10, 10, 10)  # equal time lags in minutes
  loc_err <- rep(10, 5)

  result <- move2utils:::.bm_variance(time_lag, loc_err, x, y)
  expect_type(result, "list")
  expect_named(result, c("BMvar", "cll"))
  expect_true(is.numeric(result$BMvar))
  expect_true(is.numeric(result$cll))
  expect_true(result$BMvar > 0)
})

test_that(".bm_variance with zero movement gives low variance", {
  # Stationary points — variance should be very low
  x <- c(0, 0, 0, 0, 0)
  y <- c(0, 0, 0, 0, 0)
  time_lag <- c(10, 10, 10, 10, 10)
  loc_err <- rep(10, 5)

  result <- move2utils:::.bm_variance(time_lag, loc_err, x, y)
  expect_true(result$BMvar < 1)
})

test_that(".bm_variance increases with more tortuous movement", {
  time_lag <- c(10, 10, 10, 10, 10)
  loc_err <- rep(10, 5)

  # Straight movement — leave-one-out deviation is low
  r_straight <- move2utils:::.bm_variance(time_lag, loc_err,
                                       x = c(0, 100, 200, 300, 400),
                                       y = c(0, 0, 0, 0, 0))
  # Tortuous movement — leave-one-out deviation is high
  r_tortuous <- move2utils:::.bm_variance(time_lag, loc_err,
                                       x = c(0, 100, 0, 100, 0),
                                       y = c(0, 100, 0, 100, 0))
  expect_true(r_tortuous$BMvar > r_straight$BMvar)
})

test_that(".bm_variance errors on unequal input lengths", {
  expect_error(
    move2utils:::.bm_variance(c(10, 10), rep(10, 5), 1:5, 1:5),
    "Not an equal number"
  )
})

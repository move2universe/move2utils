test_that(".delta_para_orth with movement along direction gives pure parallel", {
  mu <- matrix(c(0, 0), ncol = 2)
  direction <- matrix(c(10, 0), ncol = 2)
  point <- matrix(c(5, 0), ncol = 2)

  result <- move2utils:::.delta_para_orth(mu, direction, point)
  expect_equal(unname(result[, "deltaOrth"]), 0, tolerance = 1e-10)
  expect_equal(unname(result[, "deltaPara"]), 5, tolerance = 1e-10)
})

test_that(".delta_para_orth with movement perpendicular to direction gives pure orthogonal", {
  mu <- matrix(c(0, 0), ncol = 2)
  direction <- matrix(c(10, 0), ncol = 2)
  point <- matrix(c(0, 5), ncol = 2)

  result <- move2utils:::.delta_para_orth(mu, direction, point)
  expect_equal(unname(result[, "deltaPara"]), 0, tolerance = 1e-10)
  expect_equal(unname(result[, "deltaOrth"]), 5, tolerance = 1e-10)
})

test_that(".delta_para_orth at 45 degrees splits equally", {
  mu <- matrix(c(0, 0), ncol = 2)
  direction <- matrix(c(10, 0), ncol = 2)
  point <- matrix(c(5, 5), ncol = 2)

  result <- move2utils:::.delta_para_orth(mu, direction, point)
  expect_equal(unname(result[, "deltaPara"]), unname(result[, "deltaOrth"]),
               tolerance = 1e-10)
  total <- sqrt(result[, "deltaPara"]^2 + result[, "deltaOrth"]^2)
  expect_equal(unname(total), sqrt(50), tolerance = 1e-10)
})

test_that(".delta_para_orth handles multiple points", {
  mu <- matrix(c(0, 0, 0, 0), ncol = 2, byrow = TRUE)
  direction <- matrix(c(10, 0, 10, 0), ncol = 2, byrow = TRUE)
  point <- matrix(c(5, 0, 0, 5), ncol = 2, byrow = TRUE)

  result <- move2utils:::.delta_para_orth(mu, direction, point)
  expect_equal(nrow(result), 2)
  expect_equal(unname(result[1, "deltaOrth"]), 0, tolerance = 1e-10)
  expect_equal(unname(result[2, "deltaPara"]), 0, tolerance = 1e-10)
})

test_that(".delta_para_orth handles same location (mu == point)", {
  mu <- matrix(c(5, 5), ncol = 2)
  direction <- matrix(c(10, 0), ncol = 2)
  point <- matrix(c(5, 5), ncol = 2)

  result <- move2utils:::.delta_para_orth(mu, direction, point)
  expect_equal(unname(result[, "deltaPara"]), 0)
  expect_equal(unname(result[, "deltaOrth"]), 0)
})

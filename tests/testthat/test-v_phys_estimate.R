test_that("returns a numeric scalar in m/s with attributes", {
  res <- v_phys_estimate(5, "flying")
  expect_type(as.numeric(res), "double")
  expect_length(as.numeric(res), 1L)
  expect_s3_class(res, "v_phys_estimate")
  expect_true(!is.null(attr(res, "ci")))
  expect_length(attr(res, "ci"), 2L)
  expect_true(attr(res, "ci")[1] <= as.numeric(res))
  expect_true(attr(res, "ci")[2] >= as.numeric(res))
  expect_equal(attr(res, "kmh"), as.numeric(res) * 3.6, tolerance = 1e-9)
  expect_identical(attr(res, "mode"), "flying")
  expect_equal(attr(res, "mass"), 5)
})

test_that("rejects invalid input", {
  expect_error(v_phys_estimate("five",     "flying"),
               class = "move2utils_v_phys_estimate_bad_mass")
  expect_error(v_phys_estimate(c(1, 2),    "flying"),
               class = "move2utils_v_phys_estimate_bad_mass")
  expect_error(v_phys_estimate(0,          "flying"),
               class = "move2utils_v_phys_estimate_bad_mass")
  expect_error(v_phys_estimate(-1,         "flying"),
               class = "move2utils_v_phys_estimate_bad_mass")
  expect_error(v_phys_estimate(NA_real_,   "flying"),
               class = "move2utils_v_phys_estimate_bad_mass")
  ## "swooping" is rejected by match.arg() inside the function -- not
  ## ours to class.
  expect_error(v_phys_estimate(5, "swooping"))
  expect_error(v_phys_estimate(5, "flying", ci_level = 0),
               class = "move2utils_v_phys_estimate_bad_ci_level")
  expect_error(v_phys_estimate(5, "flying", ci_level = 1),
               class = "move2utils_v_phys_estimate_bad_ci_level")
  expect_error(v_phys_estimate(5, "flying", ci_level = -0.1),
               class = "move2utils_v_phys_estimate_bad_ci_level")
})

test_that("warns on extrapolation", {
  expect_warning(v_phys_estimate(1e-10, "flying"), "outside the range")
  expect_warning(v_phys_estimate(1e6,   "running"), "outside the range")
})

test_that("Hirt 2017 example predictions match expected order of magnitude", {
  ## golden eagle ~5 kg flying: paper-level expectation 25-35 m/s
  v_eagle <- as.numeric(v_phys_estimate(5, "flying"))
  expect_gt(v_eagle, 20)
  expect_lt(v_eagle, 40)

  ## red fox ~6 kg running: paper-level expectation 13-17 m/s
  v_fox <- as.numeric(v_phys_estimate(6, "running"))
  expect_gt(v_fox, 10)
  expect_lt(v_fox, 20)

  ## bottlenose dolphin ~250 kg swimming: 8-15 m/s
  v_dolphin <- as.numeric(v_phys_estimate(250, "swimming"))
  expect_gt(v_dolphin, 5)
  expect_lt(v_dolphin, 20)
})

test_that("predictions are monotonic with mass on the small-body rising flank", {
  ## On the rising flank of the hump-shaped curve, v_max increases
  ## monotonically with mass.  The flying mode peaks at ~1 kg, so
  ## the shared-monotonicity range stops below that.
  m_seq <- c(1e-6, 1e-4, 1e-2, 0.1)
  for (mode in c("flying", "running", "swimming")) {
    v_seq <- vapply(m_seq,
                    function(m) suppressWarnings(
                      as.numeric(v_phys_estimate(m, mode))),
                    numeric(1))
    expect_true(all(diff(v_seq) > 0), info = mode)
  }
})

test_that("ci_level controls the width of the parameter interval", {
  v_95 <- v_phys_estimate(5, "flying", ci_level = 0.95)
  v_50 <- v_phys_estimate(5, "flying", ci_level = 0.50)
  w_95 <- diff(attr(v_95, "ci"))
  w_50 <- diff(attr(v_50, "ci"))
  expect_true(w_95 > w_50)
})

test_that("print method runs and shows the citation", {
  res <- v_phys_estimate(5, "flying")
  expect_output(print(res), "Hirt")
  expect_output(print(res), "v_max")
})

test_that("mt_thin_distance rejects bad inputs", {
  expect_error(mt_thin_distance(data.frame(x = 1), distance = 100),
               "inherits")
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M1", ][1:10, ]
  expect_error(mt_thin_distance(x, distance = -1))
  expect_error(mt_thin_distance(x, distance = 0))
})

test_that("mt_thin_distance adds a thin_selected column of correct length", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M4", ][1:300, ]

  out <- mt_thin_distance(x, distance = 50)
  expect_s3_class(out, "move2")
  expect_true("thin_selected" %in% names(out))
  expect_type(out$thin_selected, "logical")
  expect_equal(length(out$thin_selected), nrow(x))
  expect_equal(nrow(out), nrow(x))
  expect_true(out$thin_selected[1])     # first fix always retained
})

test_that("retention count scales inversely with `distance`", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M4", ][1:300, ]
  x <- sf::st_transform(x, move2::mt_aeqd_crs(x))  # planar metres
  total_len <- sum(as.numeric(move2::mt_distance(x)), na.rm = TRUE)

  kept <- integer(0)
  for (step in c(100, 300, 1000)) {
    out <- mt_thin_distance(x, distance = step, remove = TRUE)
    ## the minimum-spacing guarantee caps retention at one fix per
    ## `step` of travel (plus the always-kept first fix).
    expect_lte(nrow(out), ceiling(total_len / step) + 1)
    expect_gte(nrow(out), 2L)
    kept <- c(kept, nrow(out))
  }
  ## larger steps retain strictly fewer fixes
  expect_true(all(diff(kept) < 0))
})

test_that("remove = TRUE returns only retained rows", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M4", ][1:300, ]

  full    <- mt_thin_distance(x, distance = 50)
  trimmed <- mt_thin_distance(x, distance = 50, remove = TRUE)
  expect_equal(nrow(trimmed), sum(full$thin_selected))
})

test_that("multi-track input processes each individual independently", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  subset_id <- function(id) fishers[move2::mt_track_id(fishers) == id, ][1:100, ]
  two <- rbind(subset_id("F1"), subset_id("F2"))

  out <- mt_thin_distance(two, distance = 50)
  by_id <- tapply(out$thin_selected,
                  droplevels(factor(move2::mt_track_id(out))),
                  sum)
  expect_true(all(by_id > 0))
  ## each track's first fix should be selected
  first_of_each <- tapply(out$thin_selected,
                          droplevels(factor(move2::mt_track_id(out))),
                          function(v) v[1])
  expect_true(all(first_of_each))
})

test_that("units input for distance is converted to metres", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M4", ][1:200, ]

  a <- mt_thin_distance(x, distance = 500)
  b <- mt_thin_distance(x, distance = units::set_units(0.5, "km"))
  expect_equal(a$thin_selected, b$thin_selected)
})

test_that("method = 'step' guarantees along-track spacing >= distance", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M1", ][1:300, ]

  sel <- mt_thin_distance(x, distance = 300, method = "step")$thin_selected
  d <- as.numeric(move2::mt_distance(x, units = "m")); d[is.na(d)] <- 0
  cumd <- cumsum(c(0, utils::head(d, -1L)))
  gaps <- diff(cumd[which(sel)])
  ## every retained pair is at least `distance` of travel apart
  expect_true(all(gaps >= 300 - 1e-6))
  expect_true(sel[1])   # first fix retained
})

test_that("method = 'interval' keeps pairs inside the tolerance band", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M1", ][1:300, ]

  D <- 300; tol <- 30
  sel <- mt_thin_distance(x, distance = D, tolerance = tol,
                          method = "interval")$thin_selected
  d <- as.numeric(move2::mt_distance(x, units = "m")); d[is.na(d)] <- 0
  cumd <- cumsum(c(0, utils::head(d, -1L)))
  gaps <- diff(cumd[which(sel)])
  ## within-band gaps lie in [D-tol, D+tol]; any larger gap is an
  ## un-bridgeable single step the run-splitter refused to cross.
  within <- gaps[gaps <= D + tol + 1e-6]
  expect_true(all(within >= D - tol - 1e-6))
  ## no retained pair is closer than D - tol
  expect_false(any(gaps < D - tol - 1e-6))
})

test_that("step is the default and differs from interval", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M1", ][1:300, ]

  default_sel  <- mt_thin_distance(x, distance = 300)$thin_selected
  step_sel     <- mt_thin_distance(x, distance = 300, method = "step")$thin_selected
  interval_sel <- mt_thin_distance(x, distance = 300, method = "interval")$thin_selected
  expect_identical(default_sel, step_sel)
  expect_false(identical(step_sel, interval_sel))
})

test_that("tolerance set under method = 'step' warns", {
  skip_if_not_installed("move2")
  fishers <- move2::mt_read(move2::mt_example())
  fishers <- fishers[!sf::st_is_empty(fishers), ]
  x <- fishers[move2::mt_track_id(fishers) == "M1", ][1:50, ]
  expect_warning(mt_thin_distance(x, distance = 300, tolerance = 30),
                 class = "move2utils_mt_thin_distance_tolerance_ignored")
})

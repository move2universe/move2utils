## Numerical-equivalence guard for the C port of .bgb_var_break.
##
## The historical pure-R implementation is preserved as
## .bgb_var_break_r_reference() so we can assert that the Brent-method
## C kernel and the optim()-based R reference converge to the same
## per-axis sigmas on real-shape windows.  Tolerance is set at the
## optimiser-convergence floor (~1e-3 absolute on sigma) -- both
## solvers find the same minimum modulo their internal tolerances.

test_that(".bgb_var_break (C) matches .bgb_var_break_r_reference on CPF windows", {
  skip_if_not_installed("move2")

  m <- move2::mt_read(system.file("extdata/synthetic_tracks.csv.gz",
                                    package = "move2utils"))
  m <- m[!sf::st_is_empty(m), ]
  m <- dplyr::arrange(m, move2::mt_track_id(m), move2::mt_time(m))
  trk <- m[move2::mt_track_id(m) == "CPF_A", ]
  trk <- sf::st_transform(trk, move2::mt_aeqd_crs(trk))

  cc <- sf::st_coordinates(trk)
  ts_min <- as.numeric(move2::mt_time(trk), units = "mins")
  loc_err <- rep(25, nrow(trk))

  ## Sample three windows: one near the start, one in the middle, one
  ## near the end.  Choose window_size = 31, margin = 5 (the same
  ## defaults used in the package's audit + sweep).
  ws <- 31L
  starts <- c(1L, floor(nrow(trk) / 2) - 15L, nrow(trk) - ws + 1L - 5L)

  for (s in starts) {
    idx <- s:(s + ws - 1L)
    ref <- move2utils:::.bgb_var_break_r_reference(
      cc[idx, 1], cc[idx, 2], ts_min[idx], loc_err[idx], margin = 5L)
    new <- move2utils:::.bgb_var_break(
      cc[idx, 1], cc[idx, 2], ts_min[idx], loc_err[idx], margin = 5L)

    expect_equal(is.na(ref$paraSd), is.na(new$paraSd))
    expect_equal(is.na(ref$orthSd), is.na(new$orthSd))
    expect_equal(ref$paraSd, new$paraSd, tolerance = 1e-3,
                 info = sprintf("paraSd mismatch at window starting %d", s))
    expect_equal(ref$orthSd, new$orthSd, tolerance = 1e-3,
                 info = sprintf("orthSd mismatch at window starting %d", s))
  }
})

test_that("mt_dbgb_variance C path runs at order-of-magnitude dBBMM speed", {
  skip_if_not_installed("move2")
  skip_on_cran()

  m <- move2::mt_read(system.file("extdata/synthetic_tracks.csv.gz",
                                    package = "move2utils"))
  m <- m[!sf::st_is_empty(m), ]
  m <- dplyr::arrange(m, move2::mt_track_id(m), move2::mt_time(m))
  trk <- m[move2::mt_track_id(m) == "CPF_A", ]
  trk <- sf::st_transform(trk, move2::mt_aeqd_crs(trk))

  t_dbbmm <- system.time(
    mt_dbbmm_variance(trk, location_error = 25, window_size = 31, margin = 5)
  )["elapsed"]
  t_dbgb <- system.time(
    mt_dbgb_variance(trk, location_error = 25, window_size = 31, margin = 5)
  )["elapsed"]

  ## Pre-port the ratio was ~200x; we expect <= 5x as a generous
  ## ceiling.  The dBGB cost is ~2x dBBMM in principle (two 1D fits
  ## per axis instead of one), with extra Brent calls for the
  ## breakpoint search; 5x leaves room for varying CI-runner CPU.
  ratio <- t_dbgb / t_dbbmm
  message(sprintf("dBGB / dBBMM ratio = %.2f (dBGB %.2fs, dBBMM %.2fs)",
                  ratio, t_dbgb, t_dbbmm))
  expect_lt(ratio, 5)
})

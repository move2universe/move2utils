## Tests for nested `pool_by` API (length-2: c(outer, inner)).
## The orchestrator's pool_by tests in test-mt_clean_track.R cover
## validation error wording.  This file focuses on the actual
## behavioural contract of the two-element form.

## ---- shared fixture --------------------------------------------------

## Build a 4-track move2 with two populations, each with two
## individuals, each as a single track.  Add a planted outlier on one
## track per population so the pool-fit behaviour has signal to detect.
make_nested_fixture <- function() {
  set.seed(123)
  build <- function(id, lon0, lat0, n = 120) {
    data.frame(id = id,
                timestamp = as.POSIXct("2026-01-01", tz = "UTC") +
                              seq_len(n) * 3600,
                lon = lon0 + cumsum(rnorm(n, 0, 0.002)),
                lat = lat0 + cumsum(rnorm(n, 0, 0.002)))
  }
  d <- rbind(build("t1", 10,   50),
              build("t2", 10.1, 50),
              build("t3", 11,   51),
              build("t4", 11.1, 51))
  ## Plant a jump-outlier in t1 and t3 (one per population) so the
  ## pool-fit cascade has something to flag.
  d$lon[d$id == "t1"][60] <- d$lon[d$id == "t1"][60] + 0.5  # ~50 km jump
  d$lon[d$id == "t3"][60] <- d$lon[d$id == "t3"][60] + 0.5
  m <- move2::mt_as_move2(d, coords = c("lon","lat"),
                            time_column = "timestamp",
                            track_id_column = "id", crs = 4326)
  td <- move2::mt_track_data(m)
  td$indv <- c("I1","I2","I3","I4")
  td$pop  <- c("P1","P1","P2","P2")
  move2::mt_set_track_data(m, td)
}


## ---- length-1 backward compatibility ---------------------------------

test_that("length-1 pool_by still works after the refactor", {
  m <- make_nested_fixture()
  expect_error(suppressMessages(suppressWarnings(
    mt_clean_track(m, pool_by = "indv",
                    plot = FALSE, remove = FALSE, silent = TRUE))),
    NA)
})

test_that("pool_by = c(x, x) is rejected with a clear message", {
  m <- make_nested_fixture()
  expect_error(
    suppressMessages(suppressWarnings(
      mt_clean_track(m, pool_by = c("indv","indv"),
                      plot = FALSE, remove = FALSE, silent = TRUE))),
    "two \\*distinct\\* columns")
})


## ---- nested form actually runs through ------------------------------

test_that("nested pool_by = c(outer, inner) runs and respects nesting", {
  m <- make_nested_fixture()
  o <- suppressMessages(suppressWarnings(
          mt_clean_track(m, pool_by = c("pop","indv"),
                          plot = FALSE, remove = FALSE, silent = TRUE)))
  expect_true("is_outlier" %in% names(o))
  ## The planted outliers should be flagged.  Locate them by track +
  ## within-track row index.
  ids <- as.character(move2::mt_track_id(o))
  for (tid in c("t1", "t3")) {
    rel <- which(ids == tid)
    expect_true(o$is_outlier[rel[60]],
                info = sprintf("planted outlier in track %s should flag",
                                tid))
  }
})


## ---- nesting violation error ---------------------------------------

test_that("nesting violation (inner spans multiple outer) errors", {
  m <- make_nested_fixture()
  ## indv nests in pop, but pop does NOT nest in indv (P1 spans I1,I2).
  expect_error(
    suppressMessages(suppressWarnings(
      mt_clean_track(m, pool_by = c("indv","pop"),
                      plot = FALSE, remove = FALSE, silent = TRUE))),
    "requires inner to nest in outer")
})


## ---- additivity preserved under nested form ------------------------

test_that("nested pool union is additive (per-track flags subset of pooled)", {
  m <- make_nested_fixture()
  o_pt <- suppressMessages(suppressWarnings(
            mt_clean_track(m, plot = FALSE, remove = FALSE, silent = TRUE)))
  o_pg <- suppressMessages(suppressWarnings(
            mt_clean_track(m, pool_by = c("pop","indv"),
                            plot = FALSE, remove = FALSE, silent = TRUE)))
  expect_true(all(which(o_pt$is_outlier) %in% which(o_pg$is_outlier)))
})


## ---- primitive-level nested wiring (bridge / detour / speed_cap) ---

test_that("primitives accept nested pool_by without erroring", {
  m <- make_nested_fixture()
  ## bridge
  expect_error(suppressMessages(suppressWarnings(
    mt_flag_outliers_bridge(m, pool_by = c("pop","indv"),
                              plot = FALSE, silent = TRUE))),
    NA)
  ## detour
  expect_error(suppressMessages(suppressWarnings(
    mt_flag_outliers_detour(m, pool_by = c("pop","indv"),
                              plot = FALSE, silent = TRUE,
                              threshold_type = "auto"))),
    NA)
  ## speed cap
  expect_error(suppressMessages(suppressWarnings(
    mt_flag_speed_cap(m, pool_by = c("pop","indv"),
                       plot = FALSE, silent = TRUE,
                       threshold_type = "auto"))),
    NA)
  ## prob
  expect_error(suppressMessages(suppressWarnings(
    mt_flag_outliers(m, pool_by = c("pop","indv"),
                      plot = FALSE, remove = FALSE, silent = TRUE))),
    NA)
})

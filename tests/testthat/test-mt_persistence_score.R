## Tests for the multi-scale persistence-score annotation helper.

read_cpf <- function(track = "CPF_A") {
  d <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")),
                 stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  x <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier", crs = 4326)
  x <- x[!sf::st_is_empty(x), ]
  if (is.null(track)) return(x)
  x[move2::mt_track_id(x) == track, ]
}

## ---- input validation ----

test_that("mt_persistence_score rejects non-move2 input", {
  expect_error(mt_persistence_score(NULL), "move2")
  expect_error(mt_persistence_score(data.frame(x = 1:3)), "move2")
})

test_that("mt_persistence_score requires is_outlier column", {
  x <- read_cpf("CPF_C")
  expect_error(mt_persistence_score(x), "is_outlier")
})

test_that("mt_persistence_score requires logical is_outlier", {
  x <- read_cpf("CPF_C")
  x$is_outlier <- 0L  # integer, not logical
  expect_error(mt_persistence_score(x), "logical")
})

test_that("mt_persistence_score validates scales", {
  x <- read_cpf("CPF_C")
  x$is_outlier <- rep(FALSE, nrow(x))
  expect_error(mt_persistence_score(x, scales = integer(0)),
               "non-empty integer vector")
  expect_error(mt_persistence_score(x, scales = c(0L, 2L)),
               ">= 2")
  expect_error(mt_persistence_score(x, scales = c(1L, 2L)),
               ">= 2")  # scale 1 is implicit
  expect_error(mt_persistence_score(x, scales = c(NA_integer_, 4L)),
               "non-empty integer vector")
})

test_that("mt_persistence_score validates threshold", {
  x <- read_cpf("CPF_C")
  x$is_outlier <- rep(FALSE, nrow(x))
  expect_error(mt_persistence_score(x, threshold = -1),
               "positive scalar")
  expect_error(mt_persistence_score(x, threshold = c(1, 2)),
               "positive scalar")
  expect_error(mt_persistence_score(x, threshold = NA_real_),
               "positive scalar")
})

test_that("mt_persistence_score validates n_breaks", {
  x <- read_cpf("CPF_C")
  x$is_outlier <- rep(FALSE, nrow(x))
  expect_error(mt_persistence_score(x, n_breaks = 2L),
               "integer >= 4")
})

## ---- schema ----

test_that("mt_persistence_score adds expected columns", {
  x <- read_cpf("CPF_A")
  res <- suppressMessages(suppressWarnings(
    mt_clean_track(x, plot = FALSE, remove = FALSE)))
  ann <- suppressMessages(mt_persistence_score(res))

  ## Default scales = c(2, 4, 8) -> three per-scale columns
  expect_true("persistence_count" %in% names(ann))
  expect_true("persistence_at_scale_2" %in% names(ann))
  expect_true("persistence_at_scale_4" %in% names(ann))
  expect_true("persistence_at_scale_8" %in% names(ann))
  expect_type(ann$persistence_count, "integer")
})

test_that("mt_persistence_score honours custom scales", {
  x <- read_cpf("CPF_A")
  res <- suppressMessages(suppressWarnings(
    mt_clean_track(x, plot = FALSE, remove = FALSE)))
  ann <- suppressMessages(mt_persistence_score(res, scales = c(3L, 6L)))

  expect_true("persistence_at_scale_3" %in% names(ann))
  expect_true("persistence_at_scale_6" %in% names(ann))
  expect_false("persistence_at_scale_2" %in% names(ann))
  expect_false("persistence_at_scale_4" %in% names(ann))
})

## ---- semantics ----

test_that("mt_persistence_score does not modify is_outlier", {
  x <- read_cpf("CPF_A")
  res <- suppressMessages(suppressWarnings(
    mt_clean_track(x, plot = FALSE, remove = FALSE)))
  before <- res$is_outlier
  ann <- suppressMessages(mt_persistence_score(res))
  expect_equal(ann$is_outlier, before)
})

test_that("mt_persistence_score persistence_count is NA on non-flagged fixes", {
  x <- read_cpf("CPF_A")
  res <- suppressMessages(suppressWarnings(
    mt_clean_track(x, plot = FALSE, remove = FALSE)))
  ann <- suppressMessages(mt_persistence_score(res))
  expect_true(all(is.na(ann$persistence_count[!ann$is_outlier])))
})

test_that("mt_persistence_score persistence_count is in [1, length(scales)+1]", {
  x <- read_cpf("CPF_A")
  res <- suppressMessages(suppressWarnings(
    mt_clean_track(x, plot = FALSE, remove = FALSE)))
  ann <- suppressMessages(mt_persistence_score(res, scales = c(2L, 4L, 8L)))
  pc <- ann$persistence_count[ann$is_outlier]
  expect_true(all(pc >= 1L & pc <= 4L))
})

test_that("mt_persistence_score handles no flagged fixes gracefully", {
  x <- read_cpf("CPF_C")
  x$is_outlier <- rep(FALSE, nrow(x))
  expect_message(
    ann <- mt_persistence_score(x),
    "no flagged fixes|defined only on candidates",
    fixed = FALSE
  )
  expect_true(all(is.na(ann$persistence_count)))
})

## ---- multi-track dispatch ----

test_that("mt_persistence_score dispatches per-track on multi-track input", {
  x <- read_cpf(NULL)  # all CPF tracks
  res <- suppressMessages(suppressWarnings(
    mt_clean_track(x, plot = FALSE, remove = FALSE)))
  ann <- suppressMessages(mt_persistence_score(res))
  ## Same number of rows; columns added; is_outlier preserved
  expect_equal(nrow(ann), nrow(res))
  expect_true("persistence_count" %in% names(ann))
  expect_equal(ann$is_outlier, res$is_outlier)
})

## ---- empirical class-conditional discrimination ----
##
## Pinning the empirical finding from the prototype: persistence
## score discriminates TPs from FPs when applied to cascade output
## within the state_anomaly + consensus error classes (TPs persist
## meaningfully more than FPs at p >= 3).  This anchors the helper's
## documented usage pattern.

test_that("persistence_count is identical for lon/lat and AEQD inputs (CRS-invariance)", {
  ## Earlier versions used raw st_coordinates for step/turn computation;
  ## this made the histogram bin breaks CRS-sensitive (the data range
  ## differs slightly between Haversine-on-lon/lat and Euclidean-on-AEQD,
  ## shifting the gap threshold).  After the 2026-05-11 fix the
  ## function auto-projects lon/lat to AEQD internally, so output is
  ## deterministic regardless of input CRS.
  x_ll <- read_cpf("CPF_D")  # comes in EPSG:4326
  out_ll <- suppressMessages(suppressWarnings(
    mt_clean_track(x_ll, plot = FALSE, remove = FALSE)))
  aeqd <- move2::mt_aeqd_crs(x_ll, center = "center", units = "m")
  x_ae <- sf::st_transform(x_ll, aeqd)
  out_ae <- suppressMessages(suppressWarnings(
    mt_clean_track(x_ae, plot = FALSE, remove = FALSE)))
  expect_identical(out_ll$is_outlier, out_ae$is_outlier)
  p_ll <- suppressMessages(mt_persistence_score(out_ll))$persistence_count
  p_ae <- suppressMessages(mt_persistence_score(out_ae))$persistence_count
  flag <- out_ll$is_outlier
  expect_identical(p_ll[flag], p_ae[flag])
})

test_that("persistence score discriminates TP from FP on cascade output", {
  ## Use CPF_D (block contamination) where the cascade has both
  ## TPs and FPs in state_anomaly + consensus classes.
  x <- read_cpf("CPF_D")
  gt <- readRDS(system.file("extdata", "synthetic_ground_truth.rds",
                              package = "move2utils"))
  truth <- gt[["CPF_D"]]$index

  ## Pinned to class_aware: this exercises persistence discrimination on
  ## the canonical conjunction taxonomy (state_anomaly + consensus
  ## classes), independent of the default flag rule.
  res <- suppressMessages(suppressWarnings(
    mt_clean_track(x, plot = FALSE, remove = FALSE,
                   consensus = "class_aware")))
  ann <- suppressMessages(mt_persistence_score(res))

  flagged <- which(ann$is_outlier)
  if (length(flagged) < 5L) skip("too few cascade flags for the discrimination test")

  is_tp <- flagged %in% truth
  pc <- ann$persistence_count[flagged]
  ## TPs in CPF_D persist on average at higher scores than FPs.
  ## Mean test (one-sided): E[pc | TP] > E[pc | FP] empirically.
  if (any(is_tp) && any(!is_tp)) {
    expect_gt(mean(pc[is_tp]), mean(pc[!is_tp]))
  }
})

## Propagation regression tests.
##
## Lock in the 2026-05-25 propagation fix
## (audits/2026-05-25-parameter-propagation/findings.md): user-
## supplied entropy_threshold / gap_threshold / persistence_filter_threshold
## must reach the leaf detectors (`.entropy_threshold_lower`,
## `.gap_threshold_lower`).
##
## Approach: monkey-patch the leaf to record the `threshold` it was
## actually called with on each invocation, then verify the user-
## supplied value (or `NULL` -> leaf-formal default) reaches the leaf
## through every documented user-facing path.  This is a direct
## propagation check; an "output must differ" check is fragile because
## the entropy/gap detectors saturate (deep valleys / large gaps
## accept any threshold inside a wide range) on most realistic data.

suppressPackageStartupMessages({
  library(move2); library(sf)
})

read_synthetic <- function() {
  path <- system.file("extdata", "synthetic_tracks.csv.gz",
                       package = "move2utils")
  if (nchar(path) == 0)
    path <- "inst/extdata/synthetic_tracks.csv.gz"
  d <- read.csv(gzfile(path), stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  m <- move2::mt_as_move2(d,
    coords = c("location.long", "location.lat"),
    time_column = "timestamp",
    track_id_column = "individual.local.identifier",
    crs = 4326)
  m <- m[!sf::st_is_empty(m), ]
  m <- dplyr::arrange(m, move2::mt_track_id(m), move2::mt_time(m))
  m
}

## Spy-and-record helper.  Patches `name` (an internal `.foo` in
## move2utils) to a wrapper that captures the value of formal `arg` on
## each call into a 0-init counter list, then forwards to the original.
##
## Use with `on.exit(restore())` to keep the namespace clean between
## tests.
##
## Returns: list with `calls` (numeric vector of captured arg values
## across calls, in invocation order) and `restore` (closure).
.install_spy <- function(name, arg) {
  ns <- asNamespace("move2utils")
  orig <- get(name, envir = ns)
  recorded <- list()
  spy <- function(...) {
    call_args <- list(...)
    recorded[[length(recorded) + 1L]] <<- call_args[[arg]]
    do.call(orig, call_args)
  }
  ## Preserve the leaf's formals so callers that pass by position still
  ## bind correctly.  Spy accepts the same signature.
  formals(spy) <- formals(orig)
  body(spy) <- bquote({
    call_args <- as.list(match.call())[-1L]
    arg_val <- if (.(arg) %in% names(call_args))
                  eval(call_args[[.(arg)]], envir = parent.frame())
                else
                  formals(sys.function())[[.(arg)]]
    recorded[[length(recorded) + 1L]] <<- arg_val
    .(body(orig))
  })
  ## The substituted body refers to `recorded` via lexical scope so we
  ## bind it via environment.
  environment(spy) <- environment()
  assignInNamespace(name, spy, ns = ns)
  list(recorded = function() recorded,
        restore  = function() assignInNamespace(name, orig, ns = ns))
}


## ---------------------------------------------------------------------
## Standalone primitive propagation.
## ---------------------------------------------------------------------

test_that("mt_flag_outliers_bridge entropy path forwards user threshold to .entropy_threshold_lower", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  s <- .install_spy(".entropy_threshold_lower", "threshold")
  on.exit(s$restore())

  suppressMessages(
    mt_flag_outliers_bridge(m_A, threshold_type = "entropy",
                              threshold = 0.42,
                              plot = FALSE, silent = TRUE))

  rec <- unlist(s$recorded())
  expect_true(length(rec) >= 1L)
  expect_true(all(rec == 0.42),
               info = paste("Captured leaf thresholds:", paste(rec, collapse = ",")))
})


test_that("mt_flag_outliers_bridge entropy path with NULL threshold reaches leaf formal default", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  s <- .install_spy(".entropy_threshold_lower", "threshold")
  on.exit(s$restore())

  suppressMessages(
    mt_flag_outliers_bridge(m_A, threshold_type = "entropy",
                              plot = FALSE, silent = TRUE))

  rec <- unlist(s$recorded())
  expect_true(length(rec) >= 1L)
  ## NULL forwards to the leaf formal (single source of truth).
  expect_true(all(rec == formals(move2utils:::.entropy_threshold_lower)$threshold),
               info = paste("Captured leaf thresholds:", paste(rec, collapse = ",")))
})


test_that("mt_flag_outliers_bridge gap path forwards user threshold to .gap_threshold_lower", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  s <- .install_spy(".gap_threshold_lower", "threshold")
  on.exit(s$restore())

  suppressMessages(
    mt_flag_outliers_bridge(m_A, threshold_type = "gap",
                              threshold = 4.2,
                              plot = FALSE, silent = TRUE))

  rec <- unlist(s$recorded())
  expect_true(length(rec) >= 1L)
  expect_true(all(rec == 4.2))
})


## ---------------------------------------------------------------------
## Cascade propagation via mt_clean_track.
## ---------------------------------------------------------------------

test_that("mt_clean_track entropy_threshold knob reaches the leaf via bridge primitive", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  s <- .install_spy(".entropy_threshold_lower", "threshold")
  on.exit(s$restore())

  suppressMessages(
    mt_clean_track(m_A, bridge_threshold_type = "entropy",
                    entropy_threshold = 0.42,
                    plot = FALSE, remove = FALSE, silent = TRUE))

  rec <- unlist(s$recorded())
  expect_true(length(rec) >= 1L)
  expect_true(0.42 %in% rec,
               info = paste("Captured leaf thresholds across calls:",
                             paste(rec, collapse = ",")))
})


test_that("mt_clean_track gap_threshold knob reaches the leaf via bridge primitive", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  s <- .install_spy(".gap_threshold_lower", "threshold")
  on.exit(s$restore())

  suppressMessages(
    mt_clean_track(m_A, bridge_threshold_type = "gap",
                    gap_threshold = 4.2,
                    plot = FALSE, remove = FALSE, silent = TRUE))

  rec <- unlist(s$recorded())
  expect_true(length(rec) >= 1L)
  expect_true(4.2 %in% rec,
               info = paste("Captured leaf thresholds:",
                             paste(rec, collapse = ",")))
})


test_that("mt_clean_track gap_threshold knob reaches the leaf via prob primitive (default prob_threshold_type = gap)", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  s <- .install_spy(".gap_threshold_lower", "threshold")
  on.exit(s$restore())

  suppressMessages(
    mt_clean_track(m_A, gap_threshold = 4.2,
                    plot = FALSE, remove = FALSE, silent = TRUE))

  rec <- unlist(s$recorded())
  expect_true(length(rec) >= 1L)
  expect_true(4.2 %in% rec,
               info = paste("Captured leaf thresholds:",
                             paste(rec, collapse = ",")))
})


test_that("mt_clean_track NULL knobs reach the leaf formal default (single source of truth)", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  s <- .install_spy(".gap_threshold_lower", "threshold")
  on.exit(s$restore())

  suppressMessages(
    mt_clean_track(m_A, plot = FALSE, remove = FALSE, silent = TRUE))

  rec <- unlist(s$recorded())
  leaf_default <- formals(move2utils:::.gap_threshold_lower)$threshold
  expect_true(length(rec) >= 1L)
  expect_true(all(rec == leaf_default),
               info = paste("Captured leaf thresholds:",
                             paste(rec, collapse = ",")))
})


## ---------------------------------------------------------------------
## assignInNamespace reachability — sweep-harness contract.
## ---------------------------------------------------------------------

test_that("assignInNamespace on .entropy_threshold_lower reaches mt_clean_track cascade", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  s <- .install_spy(".entropy_threshold_lower", "threshold")
  on.exit(s$restore())

  ## Sweep harness pattern: override the formal of the spy (which is
  ## now the namespace entry).  Subsequent calls without explicit
  ## threshold see the perturbed formal.
  ns <- asNamespace("move2utils")
  fn <- get(".entropy_threshold_lower", envir = ns)
  formals(fn)$threshold <- 0.42   # sweep override
  assignInNamespace(".entropy_threshold_lower", fn, ns = ns)

  suppressMessages(
    mt_clean_track(m_A, bridge_threshold_type = "entropy",
                    plot = FALSE, remove = FALSE, silent = TRUE))

  rec <- unlist(s$recorded())
  expect_true(length(rec) >= 1L)
  ## The cascade leaves entropy_threshold = NULL by default, so the
  ## leaf formal (now 0.42 via assignInNamespace) is the value that
  ## reaches it.
  expect_true(all(rec == 0.42),
               info = paste("Captured leaf thresholds:",
                             paste(rec, collapse = ",")))
})


test_that("assignInNamespace on .gap_threshold_lower reaches mt_clean_track cascade", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  s <- .install_spy(".gap_threshold_lower", "threshold")
  on.exit(s$restore())

  ns <- asNamespace("move2utils")
  fn <- get(".gap_threshold_lower", envir = ns)
  formals(fn)$threshold <- 4.2
  assignInNamespace(".gap_threshold_lower", fn, ns = ns)

  suppressMessages(
    mt_clean_track(m_A, plot = FALSE, remove = FALSE, silent = TRUE))

  rec <- unlist(s$recorded())
  expect_true(length(rec) >= 1L)
  expect_true(all(rec == 4.2),
               info = paste("Captured leaf thresholds:",
                             paste(rec, collapse = ",")))
})


## ---------------------------------------------------------------------
## persistence_filter_threshold (class-aware filter inside cascade).
## ---------------------------------------------------------------------

test_that("mt_clean_track persistence_filter_threshold reaches the persistence leaf", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  s <- .install_spy(".gap_threshold_lower", "threshold")
  on.exit(s$restore())

  ## persistence_filter must be class_aware for the filter pass to run.
  suppressMessages(
    mt_clean_track(m_A, persistence_filter = "class_aware",
                    persistence_filter_threshold = 4.2,
                    plot = FALSE, remove = FALSE, silent = TRUE))

  rec <- unlist(s$recorded())
  ## The cascade runs .gap_threshold_lower in several places; one of
  ## them is .flag_at_scale inside the persistence filter (called once
  ## per scale).  Confirm the persistence-supplied value appears in
  ## the captured stream.
  if (length(rec) > 0L) {
    expect_true(4.2 %in% rec,
                 info = paste("Captured leaf thresholds:",
                               paste(rec, collapse = ",")))
  } else {
    skip("persistence filter did not fire on CPF_A this run.")
  }
})


## ---------------------------------------------------------------------
## mt_flag_outliers_dbgb residual_max auto path -> leaf.
## ---------------------------------------------------------------------

test_that("mt_flag_outliers_dbgb residual_entropy_threshold reaches .entropy_threshold_lower", {
  m <- read_synthetic()
  m_A <- m[move2::mt_track_id(m) == "CPF_A", ]
  s <- .install_spy(".entropy_threshold_lower", "threshold")
  on.exit(s$restore())

  ## Trigger the auto residual_max path by leaving residual_max = NULL.
  suppressMessages(suppressWarnings(
    mt_flag_outliers_dbgb(m_A,
                            residual_entropy_threshold = 0.42,
                            plot = FALSE, silent = TRUE)))

  rec <- unlist(s$recorded())
  expect_true(length(rec) >= 1L)
  expect_true(0.42 %in% rec,
               info = paste("Captured leaf thresholds:",
                             paste(rec, collapse = ",")))
})

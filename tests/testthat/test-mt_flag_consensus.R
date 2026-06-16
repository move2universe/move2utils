## Tests for mt_flag_consensus(): externalised voting/consensus rule.

suppressPackageStartupMessages({
  library(move2); library(sf)
})

read_synthetic <- function() {
  d <- read.csv(gzfile(system.file("extdata", "synthetic_tracks.csv.gz",
                                     package = "move2utils")),
                stringsAsFactors = FALSE)
  d$timestamp <- as.POSIXct(d$timestamp, tz = "UTC")
  m <- mt_as_move2(d, coords = c("location.long", "location.lat"),
                    time_column = "timestamp",
                    track_id_column = "individual.local.identifier",
                    crs = 4326)
  m <- m[!sf::st_is_empty(m), ]
  m <- dplyr::arrange(m, mt_track_id(m), mt_time(m))
  m
}

## Build a small move2 with hand-set per-detector flags so each
## consensus mode can be checked against known expected outputs.
fake_flagged_track <- function() {
  m <- read_synthetic()
  m_C <- m[move2::mt_track_id(m) == "CPF_C", ][1:10, ]
  m_C$flagged_by_bridge <- c(FALSE, TRUE,  FALSE, TRUE,  TRUE,  TRUE,  FALSE, TRUE,  FALSE, FALSE)
  m_C$flagged_by_prob   <- c(FALSE, FALSE, TRUE,  TRUE,  FALSE, TRUE,  TRUE,  FALSE, FALSE, FALSE)
  m_C$flagged_by_speed  <- c(FALSE, FALSE, FALSE, TRUE,  FALSE, TRUE,  TRUE,  FALSE, TRUE,  FALSE)
  m_C$flagged_by_detour <- c(FALSE, TRUE,  FALSE, TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE)
  ## indices:                  1     2      3      4      5      6      7      8      9     10
  ## votes/idx:                0     2      1      4      2      3      2      1      1      0
  m_C
}

test_that("class_aware flags only on real class rules", {
  ## Pinned to class_aware (the default flipped to evidence_corroborated
  ## 2026-06-07); this guards the conjunction class-rule behaviour.
  x <- fake_flagged_track()
  out <- mt_flag_consensus(x, mode = "class_aware")
  ## class_aware rules:
  ##   consensus (votes>=3): idx 4 (4 detectors), idx 6 (3 detectors)
  ##   geometric_spike (bridge & detour): idx 2, 4
  ##   state_anomaly ((bridge|detour) & speed): idx 4, 6
  ##   kinematic_confluence ((bridge|detour) & prob): idx 2 (wait: detour=T, prob=T -> kin), 4, 6
  ## Union: idx 2, 4, 6 fire.  Idx 5 (bridge+detour, no kin) -> geometric_spike. So 2,4,5,6.
  ## Idx 3 (only prob) -> single detector, no fire.
  ## Idx 7 (prob+speed) -> no geometric -> no class fires.
  ## Idx 8 (only bridge), 9 (only speed) -> single, no fire.
  expect_equal(which(out$is_outlier), c(2, 4, 5, 6))
})

test_that("strict mode requires one geometric AND one kinematic", {
  x <- fake_flagged_track()
  out <- mt_flag_consensus(x, mode = "strict")
  ## strict: (bridge|detour) & (prob|speed)
  ## idx 1: no -> no; 2: (T) & (F) = no; 3: F&T = no; 4: T&T = yes;
  ## 5: T&F = no; 6: T&T = yes; 7: F&T = no; 8: T&F = no; 9: F&T = no; 10: no.
  expect_equal(which(out$is_outlier), c(4, 6))
})

test_that("majority flags whenever >= 2 detectors agree", {
  x <- fake_flagged_track()
  out <- mt_flag_consensus(x, mode = "majority")
  ## votes >= 2: idx 2(2), 4(4), 5(2), 6(3), 7(2)
  expect_equal(which(out$is_outlier), c(2, 4, 5, 6, 7))
})

test_that("speed_trusted flags on speed alone OR (bridge|detour)&prob", {
  x <- fake_flagged_track()
  out <- mt_flag_consensus(x, mode = "speed_trusted")
  ## speed_trusted: speed | (geo & prob)
  ## speed alone: idx 4, 6, 7, 9
  ## geo & prob: idx 2 (detour=T, prob=F=no, ah wait: prob is F at idx 2 ... let me recheck):
  ##   idx 2: bridge=T detour=T -> geo=T; prob=F -> no.
  ##   idx 3: bridge=F detour=F -> geo=F; prob=T -> no.
  ##   idx 4: speed=T -> yes (speed)
  ##   idx 6: speed=T -> yes (speed); also geo=T & prob=T -> yes
  ##   idx 7: speed=T -> yes (speed); geo=F so geo&prob=no
  ##   idx 9: speed=T -> yes (speed)
  ## So: 4, 6, 7, 9
  expect_equal(which(out$is_outlier), c(4, 6, 7, 9))
})

test_that("any flags whenever any single detector fires", {
  x <- fake_flagged_track()
  out <- mt_flag_consensus(x, mode = "any")
  ## any fix with at least one TRUE
  expect_equal(which(out$is_outlier), c(2, 3, 4, 5, 6, 7, 8, 9))
})

test_that("custom mode accepts a user function", {
  x <- fake_flagged_track()
  ## custom: only flag when bridge fires AND speed fires
  out <- mt_flag_consensus(x, mode = "custom",
                            custom = function(b, p, s, d) b & s)
  expect_equal(which(out$is_outlier), c(4, 6))
})

test_that("custom mode errors helpfully when function is missing", {
  x <- fake_flagged_track()
  expect_error(mt_flag_consensus(x, mode = "custom"),
                "requires `custom`")
})

test_that("custom function returning wrong-length vector errors", {
  x <- fake_flagged_track()
  expect_error(mt_flag_consensus(x, mode = "custom",
                                  custom = function(b, p, s, d) c(TRUE, FALSE)),
                "same length")
})

test_that("missing detector columns are treated as silent (FALSE)", {
  x <- fake_flagged_track()
  ## Drop the detour column; class_aware rules should fall back to
  ## three-detector form.
  x_no_d <- x; x_no_d$flagged_by_detour <- NULL
  out <- mt_flag_consensus(x_no_d, mode = "any")
  ## any with three detectors: idx 2(b), 3(p), 4(b+p+s), 5(b), 6(b+p+s),
  ## 7(p+s), 8(b), 9(s) -- detour at idx 5 is no longer there but bridge
  ## already covers it.  Same union mod detour-only cells.
  expect_true(out$is_outlier[2])
  expect_true(out$is_outlier[3])
  expect_true(out$is_outlier[4])
})

test_that("input validation: non-move2 errors", {
  expect_error(mt_flag_consensus(data.frame(flagged_by_bridge = TRUE)),
                "must be a move2 object")
})

test_that("custom column names are honoured", {
  x <- fake_flagged_track()
  names(x)[names(x) == "flagged_by_bridge"] <- "my_bridge_flag"
  out <- mt_flag_consensus(x, mode = "any", bridge_col = "my_bridge_flag")
  ## Should still flag rows where the renamed bridge column is TRUE.
  expect_true(out$is_outlier[2])  # bridge fires at idx 2
})

test_that("NA values in detector columns are treated as FALSE", {
  x <- fake_flagged_track()
  x$flagged_by_bridge[3] <- NA
  out <- mt_flag_consensus(x, mode = "any")
  ## idx 3 had only bridge=NA (treated FALSE) and prob=TRUE -> still flagged.
  expect_true(out$is_outlier[3])
  ## idx 1 (no flags, no NA) -> not flagged
  expect_false(out$is_outlier[1])
})

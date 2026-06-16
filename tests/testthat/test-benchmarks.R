# ---- detection accuracy tests ----

test_that("mt_flag_outliers detects injected outliers (histogram method)", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")

  crw <- .simulate_crw_move2(n = 500, step_sd = 0.001, seed = 42)
  inj <- .inject_outliers(crw, n_outliers = 8, displacement = 0.05, seed = 99)

  result <- mt_flag_outliers(inj$data, threshold = 0.01, plot = FALSE,
                             method = "histogram", threshold_type = "significance")

  flagged_idx <- which(result$is_outlier)
  ## an outlier at index i affects steps i-1 and i, so check +/- 1 neighborhood
  nearby_injected <- unique(c(inj$injected_idx,
                              inj$injected_idx - 1,
                              inj$injected_idx + 1))
  tp <- sum(flagged_idx %in% nearby_injected)

  tpr <- tp / length(flagged_idx)  ## precision: flagged that are near injected
  recall <- sum(sapply(inj$injected_idx, function(idx)
    any(flagged_idx %in% (idx + -1:1)))) / length(inj$injected_idx)

  message(sprintf("Histogram: precision = %.2f, recall = %.2f", tpr, recall))

  ## at least some injected outliers should be recovered (within +/- 1)
  expect_true(recall >= 0.3,
              info = paste0("Recall too low: ", round(recall, 2)))
})

test_that("mt_flag_outliers detects injected outliers (copula method)", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")
  skip_if_not_installed("circular")
  skip_if_not_installed("MASS")

  crw <- .simulate_crw_move2(n = 500, step_sd = 0.001, seed = 42)
  inj <- .inject_outliers(crw, n_outliers = 8, displacement = 0.05, seed = 99)

  result <- mt_flag_outliers(inj$data, threshold = 0.01, plot = FALSE,
                             method = "copula", threshold_type = "significance")

  flagged_idx <- which(result$is_outlier)
  recall <- sum(sapply(inj$injected_idx, function(idx)
    any(flagged_idx %in% (idx + -1:1)))) / length(inj$injected_idx)

  message(sprintf("Copula: recall = %.2f", recall))

  expect_true(recall >= 0.3,
              info = paste0("Recall too low: ", round(recall, 2)))
})

test_that("mt_flag_outliers detects injected outliers (different prob_types)", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")

  crw <- .simulate_crw_move2(n = 500, step_sd = 0.001, seed = 42)
  inj <- .inject_outliers(crw, n_outliers = 8, displacement = 0.05, seed = 99)

  for (pt in c("joint", "step_turn", "custom")) {
    result <- mt_flag_outliers(inj$data, threshold = 0.01, prob_type = pt,
                               plot = FALSE, threshold_type = "significance")
    flagged_idx <- which(result$is_outlier)
    ## check within +/- 1 neighborhood
    near_count <- sum(sapply(inj$injected_idx, function(idx)
      any(flagged_idx %in% (idx + -1:1))))

    message(sprintf("prob_type = %s: %d of %d injected detected (near)",
                    pt, near_count, length(inj$injected_idx)))
    ## at least some injected outliers should be recovered nearby
    expect_true(near_count >= 1,
                info = paste0("prob_type '", pt,
                              "': no injected outliers detected nearby"))
  }
})

test_that("autodiff_alpha affects joint probability", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")

  crw <- .simulate_crw_move2(n = 200, step_sd = 0.001, seed = 42)

  r05 <- mt_flag_outliers(crw, autodiff_alpha = 0.5, plot = FALSE)
  r10 <- mt_flag_outliers(crw, autodiff_alpha = 1.0, plot = FALSE)

  ## probabilities should differ when alpha differs
  jp05 <- r05$joint_prob[!is.na(r05$joint_prob)]
  jp10 <- r10$joint_prob[!is.na(r10$joint_prob)]

  expect_false(all(jp05 == jp10),
               info = "Joint probabilities identical for different alpha values")
})

test_that("reference=x (pooled) produces different results from per-individual", {
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]

  result_ind <- mt_flag_outliers(fishers, plot = FALSE)
  result_ref <- mt_flag_outliers(fishers, reference = fishers, plot = FALSE)

  ## both should have the same number of rows
  expect_equal(nrow(result_ind), nrow(result_ref))
  ## but joint probabilities should generally differ
  expect_false(
    all(result_ind$joint_prob == result_ref$joint_prob, na.rm = TRUE),
    info = "Pooled reference and per-individual results should differ"
  )
})


# ---- speed benchmarks ----

test_that("histogram method runs within time limit", {
  skip_on_cran()
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  leroy <- fishers[mt_track_id(fishers) == "M4", ]
  leroy <- leroy[!st_is_empty(leroy), ]

  elapsed <- system.time(
    mt_flag_outliers(leroy, plot = FALSE, method = "histogram")
  )["elapsed"]
  message(sprintf("Histogram method: %.2f seconds (%d locations)",
                  elapsed, nrow(leroy)))

  ## generous limit: should complete in under 30 seconds
  expect_true(elapsed < 30,
              info = paste0("Histogram method too slow: ", round(elapsed, 2), "s"))
})

test_that("copula method runs within time limit", {
  skip_on_cran()
  skip_if_not_installed("move2")
  skip_if_not_installed("lwgeom")
  skip_if_not_installed("circular")
  skip_if_not_installed("MASS")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  leroy <- fishers[mt_track_id(fishers) == "M4", ]
  leroy <- leroy[!st_is_empty(leroy), ]

  elapsed <- system.time(
    mt_flag_outliers(leroy, plot = FALSE, method = "copula")
  )["elapsed"]
  message(sprintf("Copula method: %.2f seconds (%d locations)",
                  elapsed, nrow(leroy)))

  ## copula should be faster than histogram; generous limit
  expect_true(elapsed < 30,
              info = paste0("Copula method too slow: ", round(elapsed, 2), "s"))
})


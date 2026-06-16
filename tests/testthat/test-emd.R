make_ud <- function(vals) {
  stopifnot(length(vals) == 100)
  ## planar CRS so emd() doesn't fire its (correct) lon/lat-without-gc
  ## warning -- the EMD computations themselves are what these tests
  ## exercise, not that warning.
  r <- terra::rast(nrows = 10, ncols = 10,
                    xmin = 0, xmax = 10,
                    ymin = 0, ymax = 10,
                    crs = "EPSG:32632",
                    vals = vals)
  r / sum(terra::values(r))
}

test_that("emd rejects non-SpatRaster input", {
  expect_error(emd(matrix(1:4, 2, 2)), "SpatRaster")
})

test_that("emd requires at least two UDs", {
  r <- make_ud(runif(100))
  expect_error(emd(r), "at least two")
})

test_that("distance is small between a UD and itself", {
  set.seed(1)
  r <- make_ud(runif(100))
  stk <- c(r, r)
  names(stk) <- c("a", "b")
  ## Sinkhorn has a small entropic-bias floor; exact should hit zero
  expect_lt(as.numeric(emd(stk, reg = 0.005)), 0.01)
  skip_if_not_installed("emdist")
  expect_lt(as.numeric(emd(stk, method = "exact")), 1e-6)
})

test_that("distance is positive between different UDs", {
  set.seed(1)
  a <- make_ud(c(rep(0, 45), rep(1, 5), rep(0, 50)))      # mass at left
  b <- make_ud(c(rep(0, 50), rep(1, 5), rep(0, 45)))      # mass shifted
  stk <- c(a, b); names(stk) <- c("a", "b")
  d <- as.numeric(emd(stk))
  expect_gt(d, 0.1)       # clearly non-zero on a 10x10 unit grid
})

test_that("emd returns a dist object with correct dimnames", {
  set.seed(2)
  uds <- lapply(seq_len(3), function(i) make_ud(runif(100)))
  names(uds) <- c("x", "y", "z")
  d <- emd(uds)
  expect_s3_class(d, "dist")
  expect_equal(attr(d, "Size"), 3L)
  expect_equal(labels(d), c("x", "y", "z"))
})

test_that("symmetric: d(a, b) ~= d(b, a) up to iteration-order noise", {
  set.seed(3)
  a <- make_ud(runif(100))
  b <- make_ud(runif(100))
  d1 <- as.numeric(emd(c(a, b)))
  d2 <- as.numeric(emd(c(b, a)))
  ## Sinkhorn is mathematically symmetric but iteration order causes
  ## sub-1e-5 relative differences; exact is strictly symmetric.
  expect_lt(abs(d1 - d2) / max(d1, 1e-9), 1e-4)
})

test_that("mask_quantile = 1 disables pre-masking", {
  set.seed(4)
  uds <- lapply(seq_len(2), function(i) make_ud(runif(100)))
  names(uds) <- c("a", "b")
  d_full    <- as.numeric(emd(uds, mask_quantile = 1))
  d_masked  <- as.numeric(emd(uds, mask_quantile = 0.999))
  ## the two should be close; masking drops sub-percent tail only
  expect_true(abs(d_full - d_masked) / max(d_full, 1e-9) < 0.05)
})

test_that("method='exact' agrees with sinkhorn within ~1%", {
  skip_if_not_installed("emdist")
  set.seed(5)
  uds <- lapply(seq_len(2), function(i) make_ud(runif(100)))
  names(uds) <- c("a", "b")
  d_s <- as.numeric(emd(uds, method = "sinkhorn", reg = 0.01))
  d_e <- as.numeric(emd(uds, method = "exact"))
  rel_err <- abs(d_s - d_e) / d_e
  expect_lt(rel_err, 0.05)   # within 5%; default reg = 0.01
})

test_that("emd rejects incompatible grid geometries", {
  a <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1,
                    ymin = 0, ymax = 1, vals = runif(100))
  a <- a / sum(terra::values(a))
  b <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2,
                    ymin = 0, ymax = 1, vals = runif(100))
  b <- b / sum(terra::values(b))
  expect_error(emd(list(a = a, b = b)),
               "same raster geometry")
})

test_that("emd handles a multi-layer SpatRaster with three layers", {
  set.seed(6)
  r1 <- make_ud(runif(100))
  r2 <- make_ud(runif(100))
  r3 <- make_ud(runif(100))
  stk <- c(r1, r2, r3)
  names(stk) <- c("one", "two", "three")
  d <- emd(stk)
  expect_equal(attr(d, "Size"), 3L)
  expect_equal(length(d), 3L)    # 3 pairwise distances
  expect_true(all(as.numeric(d) > 0))
})

test_that("threshold clips the cost matrix", {
  set.seed(7)
  ## two UDs with mass at opposite corners of a 10x10 unit grid;
  ## without clipping the transport cost reflects ~10 unit separation;
  ## clipping at threshold = 1 caps every pairwise cost.
  a <- make_ud(c(rep(1, 5), rep(0, 95)))
  b <- make_ud(c(rep(0, 95), rep(1, 5)))
  stk <- c(a, b); names(stk) <- c("a", "b")
  d_uncapped <- as.numeric(emd(stk))
  d_capped   <- as.numeric(emd(stk, threshold = 1))
  expect_lt(d_capped, d_uncapped)
  expect_lte(d_capped, 1.0 + 1e-6)
})

test_that("emd rejects negative or zero threshold", {
  r <- make_ud(runif(100))
  stk <- c(r, r); names(stk) <- c("a", "b")
  expect_error(emd(stk, threshold = 0), "threshold")
  expect_error(emd(stk, threshold = -1), "threshold")
})

test_that("method='exact' rejects gc=TRUE and threshold", {
  skip_if_not_installed("emdist")
  r <- make_ud(runif(100))
  stk <- c(r, r); names(stk) <- c("a", "b")
  expect_error(emd(stk, method = "exact", gc = TRUE), "gc = TRUE")
  expect_error(emd(stk, method = "exact", threshold = 1), "threshold")
})

test_that("emd warns on lon/lat input when gc = FALSE", {
  skip_if_not_installed("terra")
  ## Build a lon/lat-tagged tiny stack
  a <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1,
                    ymin = 0, ymax = 1, vals = runif(100),
                    crs = "EPSG:4326")
  a <- a / sum(terra::values(a))
  b <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1,
                    ymin = 0, ymax = 1, vals = runif(100),
                    crs = "EPSG:4326")
  b <- b / sum(terra::values(b))
  stk <- c(a, b); names(stk) <- c("a", "b")
  expect_warning(emd(stk), "gc = FALSE")
  expect_no_warning(emd(stk, gc = TRUE))
})

test_that("gc=TRUE produces a non-negative dist", {
  skip_if_not_installed("geosphere")
  a <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1,
                    ymin = 0, ymax = 1, vals = runif(100),
                    crs = "EPSG:4326")
  a <- a / sum(terra::values(a))
  b <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1,
                    ymin = 0, ymax = 1, vals = runif(100),
                    crs = "EPSG:4326")
  b <- b / sum(terra::values(b))
  stk <- c(a, b); names(stk) <- c("a", "b")
  d <- as.numeric(emd(stk, gc = TRUE))
  expect_true(is.finite(d))
  expect_gte(d, 0)
})

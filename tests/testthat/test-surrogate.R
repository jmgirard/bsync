library(testthat)

# =========================================================================
# --- 1. GENERATOR TESTS --------------------------------------------------
# =========================================================================

test_that("Surrogate generators return matrices of correct dimensions", {
  y <- 1:50
  n_surr <- 15

  surr_circ <- generate_surrogate_circular(y, n_surrogates = n_surr)
  expect_true(is.matrix(surr_circ))
  expect_equal(dim(surr_circ), c(50, 15))

  surr_phase <- generate_surrogate_phase(y, n_surrogates = n_surr)
  expect_true(is.matrix(surr_phase))
  expect_equal(dim(surr_phase), c(50, 15))
})

test_that("Circular shift strictly preserves data distribution", {
  y <- rnorm(50)
  surr_circ <- generate_surrogate_circular(y, n_surrogates = 5)

  expect_equal(sum(surr_circ[, 1]), sum(y))
  expect_equal(mean(surr_circ[, 3]), mean(y))
  expect_true(surr_circ[1, 1] %in% y)
})

test_that("Phase randomization preserves mean and variance", {
  y <- rnorm(100, mean = 5, sd = 2)
  surr_phase <- generate_surrogate_phase(y, n_surrogates = 5)

  expect_equal(mean(surr_phase[, 1]), mean(y), tolerance = 1e-10)
  expect_equal(var(surr_phase[, 1]), var(y), tolerance = 1e-10)
})

test_that("generate_surrogate_circular warns and errors correctly", {
  y <- 1:15
  expect_error(
    generate_surrogate_circular(y, n_surrogates = 10, lag_max = 10),
    "too short relative to"
  )

  expect_warning(
    generate_surrogate_circular(1:30, n_surrogates = 25, lag_max = 5),
    "Limited unique shifts available"
  )
})

# =========================================================================
# --- 2. INTEGRATION TESTS (REAL PIPELINE) --------------------------------
# =========================================================================

test_that("Evaluators error on bad surrogate matrix inputs", {
  x <- 1:20
  y <- 1:20

  # Loosened the regex to ignore cli backtick formatting
  expect_error(
    wcc_surrogate(
      x,
      y,
      y_surrogates = c(1, 2, 3),
      window_size = 5,
      lag_max = 2
    ),
    "must be a matrix"
  )

  expect_error(
    wdtw_surrogate(
      x,
      y,
      y_surrogates = c(1, 2, 3),
      window_size = 5,
      lag_max = 2
    ),
    "must be a matrix"
  )

  bad_mat <- matrix(0, nrow = 10, ncol = 5)
  expect_error(
    wgranger_surrogate(x, y, y_surrogates = bad_mat, window_size = 5),
    "same number of rows"
  )
})

test_that("WCC surrogate pipeline integrates and returns valid object", {
  x <- rnorm(50)
  y <- rnorm(50)
  y_surr_mat <- generate_surrogate_circular(y, n_surrogates = 5, lag_max = 5)

  res <- wcc_surrogate(
    x,
    y,
    y_surrogates = y_surr_mat,
    window_size = 10,
    lag_max = 5
  )

  expect_s3_class(res, "wcc_surr")
  expect_type(res, "list")
  expect_equal(length(res$surrogate_z), 5)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

test_that("WDTW surrogate pipeline integrates and returns valid object", {
  x <- rnorm(50)
  y <- rnorm(50)
  y_surr_mat <- generate_surrogate_circular(y, n_surrogates = 5, lag_max = 5)

  # FIX 1: lag_max = 5 added
  res <- wdtw_surrogate(
    x,
    y,
    y_surrogates = y_surr_mat,
    window_size = 10,
    lag_max = 5
  )

  expect_s3_class(res, "wdtw_surr")
  expect_type(res, "list")
  expect_equal(length(res$surrogate_cost), 5)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

test_that("WGranger surrogate pipeline integrates and returns valid object", {
  x <- rnorm(50)
  y <- rnorm(50)
  y_surr_mat <- generate_surrogate_circular(y, n_surrogates = 5, lag_max = 5)

  # FIX 2: ar_order = 1 used instead of order = 1
  res <- wgranger_surrogate(
    x,
    y,
    y_surrogates = y_surr_mat,
    window_size = 10,
    ar_order = 1
  )

  expect_s3_class(res, "wgranger_surr")
  expect_type(res, "list")
  expect_equal(length(res$surrogate_f_xy), 5)
  expect_true(res$p_value_xy >= 0 && res$p_value_xy <= 1)
})

# =========================================================================
# --- 2b. ADDITIONAL ROBUSTNESS TESTS -------------------------------------
# =========================================================================

test_that("wcc_surrogate handles NA values in surrogate matrices without crashing", {
  set.seed(7)
  x <- rnorm(50)
  y <- rnorm(50)
  y_surr_mat <- generate_surrogate_circular(y, n_surrogates = 5, lag_max = 5)

  # Inject NAs into one surrogate column
  y_surr_mat[10:15, 2] <- NA

  # na.rm = TRUE (default): should return finite p-value
  res_true <- wcc_surrogate(
    x, y, y_surrogates = y_surr_mat,
    window_size = 10, lag_max = 5, na.rm = TRUE
  )
  expect_true(res_true$p_value >= 0 && res_true$p_value <= 1)

  # na.rm = FALSE: surrogate 2 will have NA-affected windows; p-value still valid
  res_false <- wcc_surrogate(
    x, y, y_surrogates = y_surr_mat,
    window_size = 10, lag_max = 5, na.rm = FALSE
  )
  expect_true(res_false$p_value >= 0 && res_false$p_value <= 1)
})

test_that("wcc_surrogate p-value is small for strongly coupled series", {
  set.seed(123)
  n <- 100
  x <- sin(seq(0, 4 * pi, length.out = n))
  y <- x + rnorm(n, sd = 0.05)  # nearly identical to x
  y_surr_mat <- generate_surrogate_circular(y, n_surrogates = 50, lag_max = 5)

  res <- wcc_surrogate(
    x, y, y_surrogates = y_surr_mat,
    window_size = 20, lag_max = 5
  )

  # Observed synchrony is real; surrogate distribution should be much lower
  expect_true(res$observed_z > mean(res$surrogate_z))
  expect_true(res$p_value <= 0.10)
})

test_that("wcc_surrogate p-value is large for independent series", {
  set.seed(456)
  n <- 100
  x <- rnorm(n)
  y <- rnorm(n)  # independent of x
  y_surr_mat <- generate_surrogate_circular(y, n_surrogates = 50, lag_max = 5)

  res <- wcc_surrogate(
    x, y, y_surrogates = y_surr_mat,
    window_size = 20, lag_max = 5
  )

  # No real synchrony; observed z should not consistently exceed surrogates
  expect_true(res$p_value > 0.05)
})

# =========================================================================
# --- 3. S3 PRINT METHOD TESTS --------------------------------------------
# =========================================================================

test_that("Print methods return silently and output text", {
  mock_wcc_obj <- list(
    observed_z = 0.8,
    surrogate_z = c(0.1, 0.2),
    p_value = 0.01,
    n_surrogates = 100
  )
  class(mock_wcc_obj) <- c("wcc_surr", "list")

  mock_wdtw_obj <- list(
    observed_cost = 10,
    surrogate_cost = c(20, 25),
    p_value = 0.6,
    n_surrogates = 100
  )
  class(mock_wdtw_obj) <- c("wdtw_surr", "list")

  expect_invisible(print(mock_wcc_obj))

  # FIX 3: expect_message instead of capture.output
  expect_message(print(mock_wcc_obj), "significantly greater")
  expect_message(print(mock_wdtw_obj), "not significantly different")
  expect_message(print(mock_wcc_obj), "too few for stable p-values")
})

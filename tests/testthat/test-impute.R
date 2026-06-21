library(testthat)

# Create a context for the test file
context("Time Series Imputation")

test_that("Mathematical imputation is correct for small gaps", {
  # Linear interpolation
  v_linear <- c(2, NA, 6, NA, 10)
  expect_equal(
    as.numeric(impute_ts_gaps(v_linear, method = "linear")),
    c(2, 4, 6, 8, 10)
  )

  # Spline interpolation (checking curvature, not just straight lines)
  # A spline through these points should dip down and curve back up
  v_spline <- c(10, 5, NA, 5, 10)
  res_spline <- as.numeric(impute_ts_gaps(v_spline, method = "spline"))

  expect_true(res_spline[3] < 5) # The interpolated minimum should be lower than 5
  expect_equal(res_spline[1], 10)
  expect_equal(res_spline[5], 10)
})

test_that("Leading and trailing NAs are never extrapolated", {
  v_edges <- c(NA, NA, 3, NA, 5, NA, NA)
  res_edges <- as.numeric(impute_ts_gaps(v_edges, method = "linear"))

  # The inner gap should be imputed
  expect_equal(res_edges[4], 4)
  # The outer edges must remain NA
  expect_true(all(is.na(res_edges[1:2])))
  expect_true(all(is.na(res_edges[6:7])))
})

test_that("Maxgap argument strictly controls imputation", {
  v_maxgap <- c(1, NA, NA, 4, NA, NA, NA, 8)

  # Maxgap of 2 should impute the first gap but leave the second
  expect_warning(
    res_maxgap <- impute_ts_gaps(v_maxgap, method = "linear", maxgap = 2),
    "Found 1 gap\\(s\\) exceeding maxgap"
  )

  res_numeric <- as.numeric(res_maxgap)
  expect_equal(res_numeric[1:4], c(1, 2, 3, 4))
  expect_true(all(is.na(res_numeric[5:7])))
  expect_equal(res_numeric[8], 8)
})

test_that("Programmatic inputs trigger appropriate errors and warnings", {
  # Non-numeric input
  expect_error(impute_ts_gaps(c("a", "b", "c")), "must be numeric")

  # Too few valid data points
  v_sparse <- c(NA, NA, 1, NA, NA)
  expect_warning(
    impute_ts_gaps(v_sparse),
    "Not enough valid data points"
  )
})

test_that("Metadata attributes are generated correctly", {
  v_meta <- c(2, NA, 6, NA, NA, NA, 10)

  # Suppress the warning for the test to cleanly check attributes
  suppressWarnings({
    res_meta <- impute_ts_gaps(v_meta, method = "linear", maxgap = 2)
  })

  summary_attr <- attr(res_meta, "imputation_summary")

  expect_type(summary_attr, "list")
  expect_equal(summary_attr$method, "linear")
  expect_equal(summary_attr$maxgap_used, 2)
  expect_equal(summary_attr$values_imputed, 1) # Only the first NA is imputed
  expect_equal(summary_attr$values_left_na, 3) # The run of 3 NAs is left alone
})

test_that("Data with no missing values returns identical vector and correct metadata", {
  v_clean <- c(1, 2, 3, 4, 5)
  res_clean <- impute_ts_gaps(v_clean)

  expect_equal(as.numeric(res_clean), v_clean)

  summary_attr <- attr(res_clean, "imputation_summary")
  expect_equal(summary_attr$values_imputed, 0)
  expect_equal(summary_attr$values_left_na, 0)
})

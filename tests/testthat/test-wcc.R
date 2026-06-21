# Tests for wcc() core functionality --------------------------------------

test_that("r_to_z handles standard and extreme bounds", {
  # Correlation of 0 should equal a Z of 0
  expect_equal(r_to_z(0), 0)

  # The function clamps 1.0 at 0.99 to prevent infinite values
  expected_max <- 0.5 * log((1 + 0.99) / (1 - 0.99))
  expect_equal(r_to_z(1), expected_max)
  expect_equal(r_to_z(-1), -expected_max)
})

test_that("wcc returns correct structure and perfect correlation at lag 0", {
  # Create two identical sine waves
  x <- sin(seq(0, 4 * pi, length.out = 100))
  y <- x

  res <- wcc(x, y, window_size = 10, lag_max = 5)

  # Check object classes and structure
  expect_s3_class(res, "wcc_res")
  expect_type(res$fisher_z, "double")
  expect_true(is.data.frame(res$results_df))

  # At lag 0 for identical series, the correlation should be exactly 1
  lag_zero_res <- res$results_df |>
    dplyr::filter(tau == 0)

  # Use > 0.99 to account for minor floating point rounding differences
  expect_true(all(lag_zero_res$wcc > 0.99))
})

test_that("wcc handles missing values without crashing", {
  x <- runif(50)
  y <- runif(50)

  # Introduce missing values
  x[10:15] <- NA

  # Test with na.rm = TRUE
  res_na_rm <- wcc(x, y, window_size = 10, lag_max = 2, na.rm = TRUE)

  # It should successfully return a valid data frame and a numeric summary
  expect_true(is.data.frame(res_na_rm$results_df))
  expect_type(res_na_rm$fisher_z, "double")
})

test_that("wcc handles time mapping correctly", {
  x <- runif(20)
  y <- runif(20)
  t_vec <- seq(0, 1.9, by = 0.1)

  res <- wcc(x, y, time = t_vec, window_size = 5, lag_max = 2)

  # The 'i' column should now contain values from t_vec, not raw indices
  expect_true(all(res$results_df$i %in% t_vec))
  expect_equal(res$settings$has_time, TRUE)
})


# Tests for wcc() input validation ----------------------------------------

test_that("wcc catches invalid inputs", {
  x <- 1:10
  y <- 1:10

  # Type and length errors for x and y
  expect_error(
    wcc(x = c("a", "b"), y = y, window_size = 2, lag_max = 1),
    "must be a numeric vector"
  )
  expect_error(
    wcc(x = x, y = c("a", "b"), window_size = 2, lag_max = 1),
    "must be a numeric vector"
  )
  expect_error(
    wcc(x = 1:5, y = 1:10, window_size = 2, lag_max = 1),
    "same length"
  )

  # Time vector errors
  expect_error(
    wcc(x, y, time = c("a", "b"), window_size = 2, lag_max = 1),
    "must be a numeric vector"
  )
  expect_error(
    wcc(x, y, time = 1:5, window_size = 2, lag_max = 1),
    "same length as"
  )

  # Hyperparameter errors
  expect_error(
    wcc(x, y, window_size = -1, lag_max = 1),
    "single positive integer"
  )
  expect_error(
    wcc(x, y, window_size = 2, lag_max = -1),
    "single positive integer"
  )
  expect_error(
    wcc(x, y, window_size = 2, lag_max = 1, window_increment = 0),
    "single positive integer"
  )
  expect_error(
    wcc(x, y, window_size = 2, lag_max = 1, lag_increment = 0),
    "single positive integer"
  )
  expect_error(
    wcc(x, y, window_size = 2, lag_max = 1, na.rm = "TRUE"),
    "single logical value"
  )
})


# Tests for suggest_wcc_params() ------------------------------------------

test_that("suggest_wcc_params calculates values correctly and warns if lag is too long", {
  # Standard case: 30Hz, 2s duration, 1s delay
  # Suppress messages to keep the test output clean
  params <- suppressMessages(
    suggest_wcc_params(
      sample_rate = 30,
      event_duration_sec = 2,
      max_delay_sec = 1,
      overlap_pct = 0.5
    )
  )

  expect_equal(params$window_size, 240)
  expect_equal(params$lag_max, 30)
  expect_equal(params$window_increment, 120)

  # Warning case: lag exceeds half the window size
  expect_warning(
    suppressMessages(suggest_wcc_params(
      sample_rate = 30,
      event_duration_sec = 1,
      max_delay_sec = 5
    )),
    "Capping `lag_max` at half the `window_size`"
  )
})


# Tests for S3 Methods ----------------------------------------------------

test_that("print and summary methods work for wcc_res objects", {
  x <- runif(20)
  y <- runif(20)
  res <- wcc(x, y, window_size = 5, lag_max = 2)

  # cli functions emit messages, not standard output
  expect_message(print(res), "Windowed Cross-Correlation Analysis")

  # Summary method header also uses cli
  expect_message(summary(res), "Cross-Correlation Value Distribution")

  # Summary method with NAs to trigger the alert message
  # We mock the results_df to contain an NA
  res_na <- res
  res_na$results_df$wcc[1] <- NA
  expect_message(summary(res_na), "missing value")
})

test_that("wcc.cpp handles edge cases and NA triggers", {
  # Generate short series for easier manual edge-case construction
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)

  # 1. Trigger Bounds Safety Check (returning NA_REAL)
  # window_size = 10, lag_max = 5 exceeds the length of the vector (5).
  # We use lag_max = 1 to satisfy the positive integer assertion.
  res_bounds <- wcc(x, y, window_size = 10, lag_max = 1)
  expect_true(all(is.na(res_bounds$results_df$wcc)))

  # 2. Trigger "valid_n <= 1" (Insufficient data)
  # By setting NA in all but one index of the window.
  # We use lag_max = 1 to satisfy the positive integer assertion.
  x_sparse <- c(1, NA, NA, NA, NA)
  y_sparse <- c(1, NA, NA, NA, NA)
  res_sparse <- wcc(x_sparse, y_sparse, window_size = 4, lag_max = 1)
  expect_true(all(is.na(res_sparse$results_df$wcc)))

  # 3. Trigger "var_x <= 0" (Zero variance)
  # A constant value window results in var = 0, triggering the division check.
  x_const <- rep(5, 10)
  y_const <- rep(5, 10)
  res_const <- wcc(x_const, y_const, window_size = 5, lag_max = 1)
  expect_true(all(is.na(res_const$results_df$wcc)))
})

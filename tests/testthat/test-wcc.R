# Tests for wcc() core functionality --------------------------------------

test_that("r_to_z handles standard and extreme bounds", {
  # Correlation of 0 should equal a Z of 0
  expect_equal(r_to_z(0), 0)

  # The function clamps 1.0 at 0.9999 to prevent infinite values
  expected_max <- base::atanh(0.9999)

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

test_that("wcc aborts when series is too short for the chosen parameters", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)
  expect_error(wcc(x, y, window_size = 10, lag_max = 1), "too short")
})

test_that("wcc.cpp handles within-window NA and zero-variance edge cases", {
  # Trigger "valid_n <= 1": only one non-NA pair in window
  # Need >= 6 elements for window_size=4, lag_max=1 to produce n_r >= 1
  x_sparse <- c(1, NA, NA, NA, NA, NA)
  y_sparse <- c(1, NA, NA, NA, NA, NA)
  res_sparse <- wcc(x_sparse, y_sparse, window_size = 4, lag_max = 1)
  expect_true(all(is.na(res_sparse$results_df$wcc)))

  # Trigger "var_x <= 0": constant window → zero variance → NA
  x_const <- rep(5, 10)
  y_const <- rep(5, 10)
  res_const <- wcc(x_const, y_const, window_size = 5, lag_max = 1)
  expect_true(all(is.na(res_const$results_df$wcc)))
})

# M1 acceptance-criteria tests --------------------------------------------

test_that("M1: realized window length is exactly window_size samples", {
  set.seed(1)
  x <- rnorm(50)
  y <- rnorm(50)

  window_size <- 10
  lag_max <- 4

  res <- wcc(x, y, window_size = window_size, lag_max = lag_max)

  # With w_max = window_size - 1, n_r = floor((50 - 9 - 8) / 1) = 33
  expected_n_windows <- floor((50L - (window_size - 1L) - 2L * lag_max) / 1L)
  expect_equal(length(unique(res$results_df$i)), expected_n_windows)

  # Old (pre-M1) formula gave floor((50 - 10 - 8) / 1) = 32 — one fewer
  expect_equal(expected_n_windows, 33L)
})

test_that("M1: na.rm = FALSE returns NA for any window containing NA", {
  set.seed(42)
  x <- rnorm(30)
  y <- rnorm(30)
  x[10] <- NA  # NA in the middle; window_size=8 means windows i=4..10 touch it

  res_true  <- wcc(x, y, window_size = 8, lag_max = 3, na.rm = TRUE)
  res_false <- wcc(x, y, window_size = 8, lag_max = 3, na.rm = FALSE)

  # na.rm=TRUE: windows beyond position 10 (i >= 11) at tau=0 should be non-NA
  clean_true <- dplyr::filter(res_true$results_df, tau == 0, i >= 11)
  expect_false(any(is.na(clean_true$wcc)))

  # na.rm=FALSE: every window whose x-slice covers position 10 (i=4..10) is NA
  na_windows <- dplyr::filter(res_false$results_df, tau == 0, i >= 4, i <= 10)
  expect_true(all(is.na(na_windows$wcc)))

  # na.rm=FALSE: windows whose x-slice does NOT cover position 10 are non-NA
  clean_false <- dplyr::filter(res_false$results_df, tau == 0, i >= 11)
  expect_false(any(is.na(clean_false$wcc)))
})

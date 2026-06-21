# Setup some simple signals for testing
set.seed(42)
n_len <- 100
t <- seq(0, 4 * pi, length.out = n_len)
sig1 <- sin(t)
# sig2 is identical to sig1 but shifted by 5 indices
sig2 <- dplyr::lag(sig1, n = 5, default = 0)

test_that("wdtw returns expected structure and classes", {
  res <- wdtw(
    x = sig1,
    y = sig2,
    window_size = 20,
    lag_max = 10,
    window_increment = 5,
    lag_increment = 1,
    scale_method = "none"
  )

  expect_s3_class(res, "wdtw_res")
  expect_type(res$settings, "list")
  expect_type(res$mean_distance, "double")
  expect_s3_class(res$results_df, "data.frame")

  # Check expected columns
  expect_true(all(
    c("row", "col", "i", "tau", "dtw_dist") %in% names(res$results_df)
  ))

  # Check basic dimensions based on the parameters
  expected_lags <- length(seq(-10, 10, by = 1))
  expect_equal(length(unique(res$results_df$tau)), expected_lags)
})

test_that("wdtw correctly scales data globally and locally", {
  # Create a drastically scaled version of sig1
  sig1_large <- sig1 * 1000

  res_unscaled <- wdtw(
    sig1,
    sig1_large,
    window_size = 10,
    lag_max = 5,
    scale_method = "none"
  )
  res_global <- wdtw(
    sig1,
    sig1_large,
    window_size = 10,
    lag_max = 5,
    scale_method = "global"
  )
  res_local <- wdtw(
    sig1,
    sig1_large,
    window_size = 10,
    lag_max = 5,
    scale_method = "local"
  )

  # Scaled distances should be drastically smaller than unscaled distances for these vectors
  expect_true(res_global$mean_distance < res_unscaled$mean_distance)
  expect_true(res_local$mean_distance < res_unscaled$mean_distance)
})

test_that("wdtw computes distance metrics L1 and L2 correctly", {
  res_l2 <- wdtw(
    sig1,
    sig2,
    window_size = 10,
    lag_max = 5,
    scale_method = "none",
    distance_metric = "L2"
  )
  res_l1 <- wdtw(
    sig1,
    sig2,
    window_size = 10,
    lag_max = 5,
    scale_method = "none",
    distance_metric = "L1"
  )

  # Squared Euclidean (L2) and Manhattan (L1) should yield different mean distances for shifted waves
  expect_false(res_l2$mean_distance == res_l1$mean_distance)
})

test_that("wdtw logic correctly identifies perfect alignment", {
  # Comparing sig1 to itself should yield a minimum distance of 0 at lag 0
  res <- wdtw(
    x = sig1,
    y = sig1,
    window_size = 15,
    lag_max = 5,
    scale_method = "none" # Keep none to avoid floating point precision issues near 0
  )

  # Find the row with the minimum overall distance
  min_idx <- which.min(res$results_df$dtw_dist)
  best_lag <- res$results_df$tau[min_idx]
  min_dist <- res$results_df$dtw_dist[min_idx]

  expect_equal(best_lag, 0)
  expect_equal(min_dist, 0)
})

test_that("S3 methods for wdtw_res work without error", {
  res <- wdtw(sig1, sig2, window_size = 20, lag_max = 10)

  # Expect no errors when printing or plotting
  expect_no_error(capture.output(print(res)))
  expect_no_error(p <- plot(res))
  expect_s3_class(p, "ggplot")
})

test_that("pick_optima handles wdtw_res correctly across search methods", {
  res <- wdtw(sig1, sig1, window_size = 15, lag_max = 5, scale_method = "none")

  # 1. Test default global search for wdtw
  # It should automatically infer search_method = "global" and find_min = TRUE
  optima_global <- pick_optima(res)

  expect_s3_class(optima_global, "wdtw_optima")
  expect_equal(attr(optima_global, "search_method"), "global")
  expect_true(attr(optima_global, "find_min"))

  # For perfectly aligned identical signals, the absolute minimum must be at lag 0
  expect_true(all(optima_global$optimum_lag == 0, na.rm = TRUE))

  # 2. Test explicit local search with strict monotonicity
  optima_local <- pick_optima(
    res,
    search_method = "local",
    L_size = 3,
    strict_monotonic = TRUE
  )

  expect_s3_class(optima_local, "wdtw_optima")
  expect_equal(attr(optima_local, "search_method"), "local")
  expect_true(attr(optima_local, "find_min"))

  # The local valley search should also successfully identify lag 0
  expect_true(all(optima_local$optimum_lag == 0, na.rm = TRUE))
})

test_that("wdtw handles out-of-bounds cleanly by returning NA", {
  res <- wdtw(sig1, sig2, window_size = 90, lag_max = 20)

  # Because window_size + lag_max > length(sig), many edge windows will fail
  # the C++ bounds check and should cleanly return NA rather than crashing
  expect_true(any(is.na(res$results_df$dtw_dist)))
})

test_that("wdtw input assertions trigger appropriate errors", {
  # 1. Vector type and length assertions
  expect_error(
    wdtw(x = c("a", "b"), y = sig2, window_size = 10, lag_max = 5),
    "must be a numeric vector"
  )
  expect_error(
    wdtw(x = sig1, y = c("a", "b"), window_size = 10, lag_max = 5),
    "must be a numeric vector"
  )
  expect_error(
    wdtw(x = sig1, y = sig2[-1], window_size = 10, lag_max = 5),
    "must be the same length"
  )

  # 2. Time vector assertions
  expect_error(
    wdtw(sig1, sig2, time = c("t1", "t2"), window_size = 10, lag_max = 5),
    "must be a numeric vector"
  )
  expect_error(
    wdtw(sig1, sig2, time = 1:5, window_size = 10, lag_max = 5),
    "must be the same length as"
  )

  # 3. Argument matching assertions for the new string parameters
  expect_error(
    wdtw(sig1, sig2, window_size = 10, lag_max = 5, scale_method = "yes"),
    "should be one of"
  )
  expect_error(
    wdtw(sig1, sig2, window_size = 10, lag_max = 5, distance_metric = "L3"),
    "should be one of"
  )

  # 4. Integer and positivity assertions
  expect_error(
    wdtw(sig1, sig2, window_size = -5, lag_max = 5),
    "single positive integer"
  )
  expect_error(
    wdtw(sig1, sig2, window_size = 10.5, lag_max = 5),
    "single positive integer"
  )
  expect_error(
    wdtw(sig1, sig2, window_size = 10, lag_max = 0),
    "single positive integer"
  )
  expect_error(
    wdtw(sig1, sig2, window_size = 10, lag_max = 5, window_increment = -1),
    "single positive integer"
  )
  expect_error(
    wdtw(sig1, sig2, window_size = 10, lag_max = 5, lag_increment = 0),
    "single positive integer"
  )
})

test_that("wdtw correctly maps the optional time vector", {
  # Create a dummy time vector
  time_vec <- seq(0, 100, length.out = length(sig1))

  res_with_time <- wdtw(
    x = sig1,
    y = sig2,
    time = time_vec,
    window_size = 10,
    lag_max = 5
  )

  # The settings flag should be marked TRUE
  expect_true(res_with_time$settings$has_time)

  # The values in the 'i' column should be mapped exactly to values in the time_vec
  # rather than raw indices
  expect_true(all(res_with_time$results_df$i %in% time_vec))
})

test_that("calc_velocity_1d handles forward and backward methods", {
  t <- 1:5
  x <- c(0, 10, 20, 10, 0)

  # Forward looks ahead, so the last value should be NA
  v_fwd <- calc_velocity_1d(t, x, method = "forward")
  expect_equal(v_fwd, c(10, 10, -10, -10, NA))

  # Backward looks behind, so the first value should be NA
  v_bwd <- calc_velocity_1d(t, x, method = "backward")
  expect_equal(v_bwd, c(NA, 10, 10, -10, -10))
})

test_that("calc_speed_1d handles forward and backward methods", {
  t <- 1:5
  x <- c(0, 10, 20, 10, 0)

  # Speed is absolute, so all non-NA values should be positive
  s_fwd <- calc_speed_1d(t, x, method = "forward")
  expect_equal(s_fwd, c(10, 10, 10, 10, NA))

  s_bwd <- calc_speed_1d(t, x, method = "backward")
  expect_equal(s_bwd, c(NA, 10, 10, 10, 10))
})

test_that("calc_speed_2d handles forward and backward methods", {
  t <- 1:5
  x <- c(0, 3, 6, 9, 12)
  y <- c(0, 4, 8, 12, 16)

  s_fwd <- calc_speed_2d(t, x, y, method = "forward")
  expect_equal(s_fwd, c(5, 5, 5, 5, NA))

  s_bwd <- calc_speed_2d(t, x, y, method = "backward")
  expect_equal(s_bwd, c(NA, 5, 5, 5, 5))
})

test_that("calc_speed_3d handles forward and backward methods", {
  t <- 1:5
  x <- c(0, 10, 20, 30, 40)
  y <- rep(0, 5)
  z <- rep(0, 5)

  s_fwd <- calc_speed_3d(t, x, y, z, method = "forward")
  expect_equal(s_fwd, c(10, 10, 10, 10, NA))

  s_bwd <- calc_speed_3d(t, x, y, z, method = "backward")
  expect_equal(s_bwd, c(NA, 10, 10, 10, 10))
})

test_that("wdtw handles NA values correctly and safely skips windows", {
  sig1_na <- sig1
  sig1_na[15] <- NA

  sig2_na <- sig2
  sig2_na[20] <- NA

  # Test NA in the 'x' vector
  res_x_na <- wdtw(
    sig1_na,
    sig2,
    window_size = 10,
    lag_max = 5,
    scale_method = "none"
  )
  expect_true(any(is.na(res_x_na$results_df$dtw_dist)))

  # Test NA in the 'y' vector
  res_y_na <- wdtw(
    sig1,
    sig2_na,
    window_size = 10,
    lag_max = 5,
    scale_method = "none"
  )
  expect_true(any(is.na(res_y_na$results_df$dtw_dist)))

  # Test NA in both vectors (with scale_method = "global" to ensure base::scale handles it)
  res_both_na <- wdtw(
    sig1_na,
    sig2_na,
    window_size = 10,
    lag_max = 5,
    scale_method = "global"
  )
  expect_true(any(is.na(res_both_na$results_df$dtw_dist)))

  # Test NA in both vectors with local scaling
  res_local_na <- wdtw(
    sig1_na,
    sig2_na,
    window_size = 10,
    lag_max = 5,
    scale_method = "local"
  )
  expect_true(any(is.na(res_local_na$results_df$dtw_dist)))
})

library(testthat)
library(dplyr)

# -------------------------------------------------------------------------
# Fixtures: Generate standard mock data for tests
# -------------------------------------------------------------------------

set.seed(42)
n_len <- 100
time_vec <- seq(0, by = 0.5, length.out = n_len)

# 1. X strongly predicts Y (Y is lag-1 of X with minor noise)
sig_x <- rnorm(n_len)
sig_y_caused_by_x <- dplyr::lag(sig_x, n = 1, default = 0) +
  rnorm(n_len, sd = 0.1)

# 2. Y strongly predicts X (X is lag-2 of Y with minor noise)
sig_y <- rnorm(n_len)
sig_x_caused_by_y <- dplyr::lag(sig_y, n = 2, default = 0) +
  rnorm(n_len, sd = 0.1)


# -------------------------------------------------------------------------
# Structural & Input Tests
# -------------------------------------------------------------------------

test_that("wgranger aborts when series is too short for the chosen window_size", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)
  expect_error(wgranger(x, y, window_size = 20), "too short")
})

test_that("wgranger input assertions trigger appropriate errors", {
  # Vector type and length assertions
  expect_error(
    wgranger(x = c("a", "b"), y = sig_y, window_size = 20),
    "must be a numeric vector"
  )
  expect_error(
    wgranger(x = sig_x, y = c("a", "b"), window_size = 20),
    "must be a numeric vector"
  )
  expect_error(
    wgranger(x = sig_x, y = sig_y[-1], window_size = 20),
    "must be the same length"
  )

  # Time vector assertions
  expect_error(
    wgranger(sig_x, sig_y, time = c("t1", "t2"), window_size = 20),
    "must be a numeric vector"
  )
  expect_error(
    wgranger(sig_x, sig_y, time = 1:5, window_size = 20),
    "must be the same length as"
  )

  # Integer and positivity assertions
  expect_error(
    wgranger(sig_x, sig_y, window_size = -5),
    "single positive integer"
  )
  expect_error(
    wgranger(sig_x, sig_y, window_size = 10, ar_order = 0),
    "single positive integer"
  )
  expect_error(
    wgranger(sig_x, sig_y, window_size = 10, window_increment = -1),
    "single positive integer"
  )
})

test_that("wgranger returns expected structure and classes", {
  res <- wgranger(
    sig_x,
    sig_y,
    window_size = 20,
    ar_order = 1,
    window_increment = 5
  )

  expect_s3_class(res, "wgranger_res")
  expect_type(res$settings, "list")
  expect_s3_class(res$results_df, "data.frame")

  # Check expected columns
  expect_true(all(
    c("i", "f_xy", "p_xy", "f_yx", "p_yx") %in% names(res$results_df)
  ))

  # Check that row count aligns with windowing math
  expected_rows <- floor((n_len - 20) / 5)
  expect_equal(nrow(res$results_df), expected_rows)
})

test_that("wgranger handles the optional time vector correctly", {
  res_time <- wgranger(sig_x, sig_y, time = time_vec, window_size = 20)

  expect_true(res_time$settings$has_time)
  # Ensure the index column maps perfectly to provided timestamps
  expect_true(all(res_time$results_df$i %in% time_vec))
})

test_that("wgranger handles insufficient degrees of freedom safely", {
  # ar_order = 10 requires 10 lags for prediction, leaving only 5 effective
  # data points in a window of 15. The model has more parameters than points.
  # Expect a warning from C++ and NA results instead of a hard crash.

  expect_warning(
    res_df <- wgranger(sig_x, sig_y, window_size = 15, ar_order = 10),
    "too small"
  )

  # All statistical results should be NA
  expect_true(all(is.na(res_df$results_df$f_xy)))
})

test_that("wgranger handles NAs safely within windows", {
  sig_x_na <- sig_x
  sig_x_na[15] <- NA

  res_na <- wgranger(sig_x_na, sig_y, window_size = 10, ar_order = 1)

  # Windows overlapping index 15 should gracefully return NA
  # In a window size of 10 with increment 1, indices 6 through 15 will contain the NA
  expect_true(any(is.na(res_na$results_df$f_xy)))
  expect_true(!all(is.na(res_na$results_df$f_xy))) # Ensure unaffected windows still compute
})


# -------------------------------------------------------------------------
# Mathematical Correctness Tests
# -------------------------------------------------------------------------

test_that("wgranger mathematically detects strict X -> Y causality", {
  # Test the scenario where Y is just a delayed X.
  # X should be an incredible predictor of Y. Y should not predict X.

  res <- wgranger(sig_x, sig_y_caused_by_x, window_size = 30, ar_order = 1)
  df <- res$results_df

  # 1. Check p-values: x->y should be highly significant (< 0.05) everywhere
  expect_true(all(df$p_xy < 0.05, na.rm = TRUE))

  # 2. Check p-values: y->x should be mostly non-significant
  # (Allowing a tiny margin for random noise triggering false positives)
  prop_sig_yx <- mean(df$p_yx < 0.05, na.rm = TRUE)
  expect_true(prop_sig_yx < 0.1)

  # 3. Check F-Stats: X->Y effect size should be massively larger than Y->X
  expect_true(mean(df$f_xy, na.rm = TRUE) > mean(df$f_yx, na.rm = TRUE) * 10)
})

test_that("wgranger mathematically detects strict Y -> X causality", {
  # Test the scenario where X is a delayed Y.
  # Note: The delay is 2 frames, so ar_order must be at least 2 to detect it!

  res <- wgranger(sig_x_caused_by_y, sig_y, window_size = 30, ar_order = 2)
  df <- res$results_df

  # 1. Check p-values: y->x should be highly significant (< 0.05) everywhere
  expect_true(all(df$p_yx < 0.05, na.rm = TRUE))

  # 2. Check p-values: x->y should be mostly non-significant
  prop_sig_xy <- mean(df$p_xy < 0.05, na.rm = TRUE)
  expect_true(prop_sig_xy < 0.1)

  # 3. Check F-Stats: Y->X effect size should be massively larger
  expect_true(mean(df$f_yx, na.rm = TRUE) > mean(df$f_xy, na.rm = TRUE) * 10)
})


# -------------------------------------------------------------------------
# S3 Method Tests
# -------------------------------------------------------------------------

test_that("S3 methods for wgranger_res work without error", {
  res <- wgranger(sig_x, sig_y, window_size = 20, ar_order = 1)

  # Printing and summarizing should not throw errors
  expect_no_error(capture.output(print(res)))
  expect_no_error(capture.output(summary(res)))
})

test_that("plot.wgranger_res handles metrics and axis scaling correctly", {
  res_raw <- wgranger(sig_x, sig_y, window_size = 20, ar_order = 1)
  res_time <- wgranger(
    sig_x,
    sig_y,
    time = time_vec,
    window_size = 20,
    ar_order = 1
  )

  # 1. Default plot (F-statistic, no native time)
  p_f <- plot(res_raw)
  expect_s3_class(p_f, "ggplot")
  expect_equal(p_f$labels$y, "F-Statistic")
  expect_equal(p_f$labels$x, "Elapsed Time Window Index")

  # 2. Metric toggle to p-values
  p_p <- plot(res_raw, metric = "p")
  expect_equal(p_p$labels$y, "-log10(p-value)")
  # Expect an extra hline and text annotation layer for the p=0.05 cutoff
  expect_true(length(p_p$layers) > length(p_f$layers))

  # 3. Scaling toggle (time_step)
  p_scaled <- plot(res_raw, time_step = 0.5)
  expect_equal(p_scaled$labels$x, "Elapsed Time (Seconds)")

  # 4. Native time vector correctly preserved
  p_native <- plot(res_time)
  expect_equal(p_native$labels$x, "Elapsed Time")

  # 5. Smooth toggle adds a layer
  p_smooth <- plot(res_raw, smooth = TRUE)
  expect_true(any(sapply(p_smooth$layers, function(l) {
    inherits(l$geom, "GeomSmooth")
  })))
})

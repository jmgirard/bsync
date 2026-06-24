# Tests for evaluate_signal_power() ---------------------------------------

test_that("evaluate_signal_power returns expected structure for single vector", {
  skip_if_not_installed("gsignal")
  skip_if_not_installed("ggplot2")

  # Create a synthetic 30Hz signal (sine wave + noise)
  fs <- 30
  t <- seq(0, 10, by = 1 / fs)
  x <- sin(2 * pi * 1 * t) + rnorm(length(t), 0, 0.1)

  res <- evaluate_signal_power(
    x,
    sample_rate = fs,
    threshold = 0.95,
    plot = TRUE
  )

  # Check structure against the newly refactored variable names
  expect_type(res, "list")
  expect_true(all(
    c(
      "primary_cutoff_freq",
      "theoretical_min_rate",
      "recommended_downsample_factor",
      "recommended_target_rate",
      "plot"
    ) %in%
      names(res)
  ))

  # Check numeric types
  expect_true(is.numeric(res$primary_cutoff_freq))
  expect_true(is.numeric(res$theoretical_min_rate))
  expect_true(is.numeric(res$recommended_downsample_factor))
  expect_true(is.numeric(res$recommended_target_rate))
  expect_true(inherits(res$plot, "ggplot"))

  # Check new integer downsampling logic
  expect_true(res$recommended_downsample_factor %% 1 == 0) # Must be an integer
  expect_equal(
    res$recommended_target_rate,
    fs / res$recommended_downsample_factor
  )

  # Ensure dataset-level metrics are NOT present for a single vector
  expect_false("all_cutoffs" %in% names(res))
})

test_that("evaluate_signal_power handles multiple signals via data frame and list", {
  skip_if_not_installed("gsignal")
  skip_if_not_installed("ggplot2")

  fs <- 30
  t <- seq(0, 10, by = 1 / fs)

  # Test with a Data Frame
  # Include a non-numeric column to verify it correctly filters it out
  df <- data.frame(
    sig1 = sin(2 * pi * 1 * t) + rnorm(length(t), 0, 0.1),
    sig2 = sin(2 * pi * 1.5 * t) + rnorm(length(t), 0, 0.1),
    sig3 = sin(2 * pi * 0.5 * t) + rnorm(length(t), 0, 0.1),
    char_col = letters[1:length(t)]
  )

  res_df <- evaluate_signal_power(
    df,
    sample_rate = fs,
    threshold = 0.95,
    plot = TRUE
  )

  expect_true(all(
    c(
      "primary_cutoff_freq",
      "theoretical_min_rate",
      "recommended_downsample_factor",
      "recommended_target_rate",
      "all_cutoffs",
      "summary_stats",
      "plot"
    ) %in%
      names(res_df)
  ))
  expect_equal(length(res_df$all_cutoffs), 3)

  # Verify the primary cutoff is exactly the 95th percentile of all cutoffs
  expected_cutoff <- unname(stats::quantile(res_df$all_cutoffs, probs = 0.95))
  expect_equal(res_df$primary_cutoff_freq, expected_cutoff)

  # Verify integer math holds for dataframes
  expect_true(res_df$recommended_downsample_factor %% 1 == 0)
  expect_equal(
    res_df$recommended_target_rate,
    fs / res_df$recommended_downsample_factor
  )

  # Test with a List
  lst <- list(a = df$sig1, b = df$sig2)
  res_lst <- evaluate_signal_power(lst, sample_rate = fs, plot = FALSE)

  expect_true(all(
    c(
      "primary_cutoff_freq",
      "theoretical_min_rate",
      "recommended_downsample_factor",
      "recommended_target_rate",
      "all_cutoffs",
      "summary_stats"
    ) %in%
      names(res_lst)
  ))
  expect_equal(length(res_lst$all_cutoffs), 2)
})

test_that("evaluate_signal_power works without plotting", {
  skip_if_not_installed("gsignal")

  x <- rnorm(100)
  res <- evaluate_signal_power(x, sample_rate = 10, plot = FALSE)

  expect_false("plot" %in% names(res))
})

test_that("evaluate_signal_power catches invalid inputs", {
  x <- rnorm(100)

  # Structural input errors
  expect_error(
    evaluate_signal_power(c("a", "b"), 30),
    "must be a numeric vector, a list of numeric vectors, or a data frame"
  )
  expect_error(
    evaluate_signal_power(data.frame(a = c("a", "b")), 30),
    "contains no numeric columns"
  )
  expect_error(
    evaluate_signal_power(list(a = 1:10, b = "string"), 30),
    "All elements in the list"
  )

  # Parameter errors
  expect_error(
    evaluate_signal_power(x, sample_rate = -5),
    "single positive number"
  )
  expect_error(
    evaluate_signal_power(x, sample_rate = 30, threshold = 1.5),
    "between 0 and 1"
  )
  expect_error(
    evaluate_signal_power(x, sample_rate = 30, plot = "TRUE"),
    "single logical value"
  )
})

test_that("evaluate_signal_power handles missing data gracefully", {
  skip_if_not_installed("gsignal")

  # Total failure: No valid data
  x_bad <- c(1, 2, rep(NA, 50))
  expect_error(
    evaluate_signal_power(x_bad, sample_rate = 10),
    "No signals contained enough non-missing data"
  )

  # Partial failure: Data frame with one good column and one broken column
  df_mixed <- data.frame(
    good = rnorm(100),
    bad = c(1, 2, rep(NA, 98))
  )

  # It should drop the bad column and successfully process the good one.
  # Because only 1 signal remains, it should revert to the single-vector output format.
  res_mixed <- evaluate_signal_power(df_mixed, sample_rate = 10, plot = FALSE)
  expect_false("all_cutoffs" %in% names(res_mixed))
  expect_true(is.numeric(res_mixed$primary_cutoff_freq))
})

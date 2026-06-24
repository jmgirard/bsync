# Tests for smooth_signal() -----------------------------------------------

test_that("smooth_signal moving_average works correctly", {
  x <- c(2, 4, 6, 8, 10)
  # A window of 3 centered on the middle elements
  smoothed <- smooth_signal(x, method = "moving_average", window = 3)

  expect_true(is.na(smoothed[1]))
  expect_true(is.na(smoothed[5]))
  expect_equal(smoothed[2], mean(c(2, 4, 6)))
  expect_equal(smoothed[3], mean(c(4, 6, 8)))
})

test_that("smooth_signal sgolay and butterworth execute if gsignal is installed", {
  # Gracefully skip these tests if the user doesn't have 'gsignal' installed
  skip_if_not_installed("gsignal")

  x <- sin(seq(0, 2 * pi, length.out = 100)) + rnorm(100, 0, 0.1)

  sg_smoothed <- smooth_signal(x, method = "sgolay", window = 5, sg_order = 3)
  expect_equal(length(sg_smoothed), 100)
  expect_false(any(is.na(sg_smoothed)))

  bw_smoothed <- smooth_signal(x, method = "butterworth", bw_cutoff = 0.2)
  expect_equal(length(bw_smoothed), 100)
  expect_false(any(is.na(bw_smoothed)))
})

test_that("smooth_signal catches invalid inputs", {
  # Main input validation
  expect_error(smooth_signal(c("a", "b", "c")), "must be a numeric vector")

  # Moving average validation
  expect_error(
    smooth_signal(1:10, method = "moving_average", window = -1),
    "positive integer"
  )
  expect_error(
    smooth_signal(1:10, method = "moving_average", window = 2.5),
    "positive integer"
  )

  # Savitzky-Golay validation
  expect_error(
    smooth_signal(1:10, method = "sgolay", window = 4),
    "odd integer"
  )
  expect_error(
    smooth_signal(1:10, method = "sgolay", window = 5, sg_order = 5),
    "strictly less than"
  )

  # Butterworth validation
  expect_error(
    smooth_signal(1:10, method = "butterworth", bw_cutoff = 1.5),
    "between 0 and 1"
  )
  expect_error(
    smooth_signal(1:10, method = "butterworth", bw_cutoff = c(0.1, 0.2)),
    "single numeric value"
  )
  expect_error(
    smooth_signal(1:10, method = "butterworth", bw_order = 0),
    "positive integer"
  )
})


# Tests for aggregate_by_time() -------------------------------------------

test_that("aggregate_by_time correctly downsamples using median and mean", {
  # Simulating 30Hz data (roughly 0.033s intervals)
  # Using values where mean and median will differ
  df <- data.frame(
    time = c(0.000, 0.033, 0.066, 0.100, 0.133, 0.166),
    au12 = c(2, 4, 12, 8, 10, 30),
    character_col = c("a", "b", "c", "d", "e", "f")
  )

  # Test median aggregation (default)
  agg_med <- aggregate_by_time(
    df,
    time_var = time,
    bin_width = 0.1,
    method = "median"
  )
  expect_false("character_col" %in% names(agg_med))
  expect_equal(nrow(agg_med), 2)
  expect_equal(agg_med$time, c(0.05, 0.15))
  # Median of c(2, 4, 12) is 4; Median of c(8, 10, 30) is 10
  expect_equal(agg_med$au12, c(4, 10))

  # Test mean aggregation
  agg_mean <- aggregate_by_time(
    df,
    time_var = time,
    bin_width = 0.1,
    method = "mean"
  )
  # Mean of c(2, 4, 12) is 6; Mean of c(8, 10, 30) is 16
  expect_equal(agg_mean$au12, c(6, 16))
})

test_that("aggregate_by_time handles na.rm correctly", {
  df <- data.frame(
    time = c(0.1, 0.2, 0.3),
    val = c(2, NA, 6)
  )

  # With na.rm = TRUE (default), the median of 2 and 6 is 4
  agg_true <- aggregate_by_time(
    df,
    time_var = time,
    bin_width = 0.5,
    na.rm = TRUE
  )
  expect_equal(agg_true$val, 4)

  # With na.rm = FALSE, the median of 2, NA, 6 is NA
  agg_false <- aggregate_by_time(
    df,
    time_var = time,
    bin_width = 0.5,
    na.rm = FALSE
  )
  expect_true(is.na(agg_false$val))
})

test_that("aggregate_by_time catches invalid inputs", {
  df <- data.frame(time = 1:10, val = 1:10)

  # Data structure and type errors
  expect_error(
    aggregate_by_time(list(a = 1), time_var = time, bin_width = 1),
    "must be a data frame"
  )
  expect_error(
    aggregate_by_time(df, time_var = time, bin_width = -1),
    "single positive number"
  )
  expect_error(
    aggregate_by_time(df, time_var = time, bin_width = c(1, 2)),
    "single positive number"
  )
  expect_error(
    aggregate_by_time(df, time_var = time, bin_width = 1, na.rm = "TRUE"),
    "single logical value"
  )
  expect_error(
    aggregate_by_time(df, time_var = time, bin_width = 1, method = "mode"),
    "should be one of"
  )
})


# Tests for trim_edges() --------------------------------------------------

test_that("trim_edges works on vectors, data frames, and matrices (pad_na = FALSE)", {
  # Vector trimming
  v <- 1:10
  expect_equal(trim_edges(v, trim_length = 2), 3:8)

  # Data frame trimming
  df <- data.frame(a = 1:5, b = letters[1:5])
  df_trimmed <- trim_edges(df, trim_length = 1)
  expect_equal(nrow(df_trimmed), 3)
  expect_equal(df_trimmed$a, 2:4)

  # Matrix trimming
  mat <- matrix(1:20, nrow = 10)
  mat_trimmed <- trim_edges(mat, trim_length = 3)
  expect_equal(nrow(mat_trimmed), 4)
})

test_that("trim_edges replaces edges with NAs when pad_na = TRUE", {
  # Vector padding
  v <- 1:5
  v_padded <- trim_edges(v, trim_length = 1, pad_na = TRUE)
  expect_equal(length(v_padded), 5)
  expect_equal(v_padded, c(NA, 2, 3, 4, NA))

  # Data frame padding
  df <- data.frame(a = 1:4, b = letters[1:4])
  df_padded <- trim_edges(df, trim_length = 1, pad_na = TRUE)
  expect_equal(nrow(df_padded), 4)
  # Check that top and bottom rows are completely NA
  expect_true(all(is.na(df_padded[1, ])))
  expect_true(all(is.na(df_padded[4, ])))
  # Check that inner rows remain intact
  expect_equal(df_padded$a[2:3], 2:3)

  # Matrix padding
  mat <- matrix(1:8, nrow = 4)
  mat_padded <- trim_edges(mat, trim_length = 1, pad_na = TRUE)
  expect_equal(nrow(mat_padded), 4)
  expect_true(all(is.na(mat_padded[1, ])))
  expect_true(all(is.na(mat_padded[4, ])))
  expect_equal(mat_padded[2:3, 1], c(2, 3))
})

test_that("trim_edges catches invalid inputs and extreme trims", {
  v <- 1:5
  df <- data.frame(a = 1:5)

  # Invalid argument types
  expect_error(trim_edges(v, trim_length = -1), "single positive integer")
  expect_error(trim_edges(v, trim_length = 1.5), "single positive integer")
  expect_error(
    trim_edges(v, trim_length = 1, pad_na = "TRUE"),
    "single logical value"
  )
  expect_error(
    trim_edges(list(a = 1:5), trim_length = 1),
    "vector, matrix, or data frame"
  )

  # Edge case where trim_length exceeds data dimensions
  expect_error(
    trim_edges(v, trim_length = 3),
    "too large; it would affect all elements"
  )
  expect_error(
    trim_edges(df, trim_length = 3),
    "too large; it would affect all rows"
  )
})


# Tests for downsample_signal() -------------------------------------------

test_that("downsample_signal handles median and mean aggregation correctly", {
  x <- c(2, 4, 12, 8, 10, 30)

  # Median of c(2,4,12) is 4; Median of c(8,10,30) is 10
  expect_equal(downsample_signal(x, factor = 3, method = "median"), c(4, 10))

  # Mean of c(2,4,12) is 6; Mean of c(8,10,30) is 16
  expect_equal(downsample_signal(x, factor = 3, method = "mean"), c(6, 16))
})

test_that("downsample_signal correctly truncates incomplete windows", {
  x <- c(1, 2, 3, 4, 5) # Length 5, factor 2 means the last element is dropped
  expect_equal(downsample_signal(x, factor = 2, method = "mean"), c(1.5, 3.5))
})

test_that("downsample_signal handles NAs correctly", {
  x <- c(2, NA, 10, 4, 5, 6)

  # Window 1: c(2, NA, 10). Median with na.rm=TRUE is 6
  expect_equal(
    downsample_signal(x, factor = 3, method = "median", na.rm = TRUE),
    c(6, 5)
  )
  expect_true(is.na(downsample_signal(
    x,
    factor = 3,
    method = "median",
    na.rm = FALSE
  )[1]))
})

test_that("downsample_signal catches invalid inputs", {
  expect_error(
    downsample_signal(c("a", "b"), factor = 2),
    "must be a numeric vector"
  )
  expect_error(downsample_signal(1:10, factor = 0), "greater than 1")
  expect_error(downsample_signal(1:10, factor = 1.5), "single integer")
  expect_error(
    downsample_signal(1:10, factor = 2, na.rm = "TRUE"),
    "single logical value"
  )

  # Removed the word "factor" from the regex to avoid backtick matching issues
  expect_error(
    downsample_signal(1:2, factor = 5),
    "smaller than the downsampling"
  )
})

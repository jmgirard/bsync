# Main Functions ----------------------------------------------------------

#' Smooth a Time Series Signal
#'
#' Applies a smoothing filter to a numeric vector. Smoothing is highly recommended
#' prior to calculating velocity or running windowed cross-correlation (WCC) to
#' reduce high-frequency noise and prevent spurious correlations.
#'
#' @param x A numeric vector representing the signal to be smoothed.
#' @param method A character string specifying the smoothing method:
#'   "moving_average", "sgolay" (Savitzky-Golay), or "butterworth". Default is "sgolay".
#' @param window An integer specifying the window size. Must be an odd number
#'   for "sgolay". Best practice: Calculate this based on the expected duration
#'   of your target behavior. (e.g., A 2-second behavior sampled at 5Hz = a window of 11).
#' @param sg_order An integer specifying the polynomial order for the Savitzky-Golay filter.
#'   Must be less than `window`. Best practice: Use 2 (quadratic) to extract broad structural
#'   trends for cross-correlation, or 3 (cubic) to preserve absolute peak intensities.
#'   Orders > 3 will typically overfit to high-frequency tracking noise. Default is 3.
#' @param bw_cutoff A numeric value between 0 and 1 specifying the normalized cutoff
#'   frequency for the Butterworth filter. Default is 0.1.
#' @param bw_order An integer specifying the order of the Butterworth filter. Default is 2.
#' @param lower_bound Numeric. If provided, smoothed values below this are clamped to this value.
#' @param upper_bound Numeric. If provided, smoothed values above this are clamped to this value.
#' @return A numeric vector containing the smoothed signal, of the same length as `x`.
#' @export
smooth_signal <- function(
  x,
  method = c("sgolay", "moving_average", "butterworth"),
  window = 5,
  sg_order = 3,
  bw_cutoff = 0.1,
  bw_order = 2,
  lower_bound = NULL,
  upper_bound = NULL
) {
  if (!is.numeric(x)) {
    cli::cli_abort("{.arg x} must be a numeric vector.")
  }

  method <- match.arg(method)

  if (method == "moving_average") {
    if (!rlang::is_integerish(window, n = 1) || window <= 0) {
      cli::cli_abort(
        "{.arg window} must be a positive integer for the moving average filter."
      )
    }
    # Uses a centered moving average to prevent phase shifts
    weights <- rep(1 / window, window)
    res <- as.numeric(stats::filter(x, weights, sides = 2))
  } else if (method == "sgolay") {
    if (!rlang::is_integerish(window, n = 1) || window %% 2 == 0) {
      cli::cli_abort(
        "{.arg window} must be an odd integer for the Savitzky-Golay filter."
      )
    }
    if (!rlang::is_integerish(sg_order, n = 1) || sg_order >= window) {
      cli::cli_abort(
        "{.arg sg_order} must be strictly less than the {.arg window} size."
      )
    }
    # Swapped from signal::sgolayfilt to gsignal::sgolayfilt
    res <- as.numeric(gsignal::sgolayfilt(x, p = sg_order, n = window))
  } else if (method == "butterworth") {
    if (
      !is.numeric(bw_cutoff) ||
        length(bw_cutoff) != 1 ||
        bw_cutoff <= 0 ||
        bw_cutoff >= 1
    ) {
      cli::cli_abort(
        "{.arg bw_cutoff} must be a single numeric value between 0 and 1."
      )
    }
    if (!rlang::is_integerish(bw_order, n = 1) || bw_order <= 0) {
      cli::cli_abort("{.arg bw_order} must be a positive integer.")
    }
    # Swapped from signal to gsignal
    bf <- gsignal::butter(bw_order, bw_cutoff, type = "low")
    res <- as.numeric(gsignal::filtfilt(bf, x))
  }

  # Apply optional bounds to prevent polynomial undershoot/overshoot
  if (!is.null(lower_bound)) {
    res <- pmax(res, lower_bound, na.rm = FALSE)
  }

  if (!is.null(upper_bound)) {
    res <- pmin(res, upper_bound, na.rm = FALSE)
  }

  return(res)
}

#' Aggregate Time Series Data by Time Bins
#'
#' Efficiently downsamples time series data by aggregating values within specified
#' time bins. This is a high-level, data frame-based pipeline function.
#'
#' @details
#' **When to use this function versus `downsample_signal()`:**
#'
#' * Use `aggregate_by_time()` when working with raw behavioral tracking data
#'   (e.g., OpenFace output) that may contain irregular timestamps, dropped
#'   frames, or missing rows. By binning based on the actual time variable,
#'   this function preserves the true chronological structure of the data and
#'   correctly leaves gaps where tracking was lost. It is also ideal for
#'   processing multiple numeric columns simultaneously.
#' * Use `downsample_signal()` when working with a single, continuous numeric
#'   vector that has guaranteed regular intervals and no missing frames. That
#'   function relies on matrix reshaping and vector math, making it exceptionally
#'   fast for clean, pre-processed data.
#'
#' @param data A data frame containing the time series data.
#' @param time_var The unquoted name of the column containing time values.
#' @param bin_width A numeric value specifying the width of the time bins.
#'   This should be in the same units as your time variable (e.g., 0.1 for 100ms bins).
#' @param method A character string specifying the aggregation method:
#'   "median" or "mean". Default is "median", which is highly robust to
#'   single-frame tracking glitches.
#' @param na.rm A logical indicating whether to remove missing values when
#'   calculating the aggregate. Default is `TRUE`.
#' @return A new data frame with the downsampled time series. The time variable
#'   is updated to represent the center of each bin, and all non-numeric columns
#'   are dropped.
#' @export
aggregate_by_time <- function(
  data,
  time_var,
  bin_width,
  method = c("median", "mean"),
  na.rm = TRUE
) {
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }
  if (!is.numeric(bin_width) || length(bin_width) != 1 || bin_width <= 0) {
    cli::cli_abort("{.arg bin_width} must be a single positive number.")
  }
  if (!rlang::is_logical(na.rm, n = 1)) {
    cli::cli_abort("{.arg na.rm} must be a single logical value.")
  }

  method <- match.arg(method)

  if (method == "median") {
    agg_fun <- function(x) stats::median(x, na.rm = na.rm)
  } else {
    agg_fun <- function(x) base::mean(x, na.rm = na.rm)
  }

  data |>
    dplyr::mutate(
      .bin_center = floor({{ time_var }} / bin_width) *
        bin_width +
        (bin_width / 2)
    ) |>
    dplyr::summarise(
      dplyr::across(dplyr::where(is.numeric), agg_fun),
      .by = .bin_center
    ) |>
    dplyr::mutate({{ time_var }} := .bin_center) |>
    dplyr::select(-.bin_center)
}

#' Trim Edge Effects from Data
#'
#' Removes or masks a specified number of observations from the beginning and end of a
#' vector or data frame. This is mathematically required after applying symmetric
#' rolling filters (like Savitzky-Golay) to remove boundary artifacts.
#'
#' @param x A numeric vector, matrix, or data frame.
#' @param trim_length An integer specifying the number of observations to mask
#'   from both ends. Best practice: For a Savitzky-Golay filter, this must exactly
#'   equal (window - 1) / 2.
#' @param pad_na A logical indicating whether to replace the trimmed edges with `NA`
#'   instead of dropping them. Set to `TRUE` when using inside `dplyr::mutate()`
#'   to preserve the original vector length. Default is `FALSE`.
#' @return An object of the same class as `x` with the edges removed or masked.
#' @export
trim_edges <- function(x, trim_length, pad_na = FALSE) {
  if (!rlang::is_integerish(trim_length, n = 1) || trim_length <= 0) {
    cli::cli_abort("{.arg trim_length} must be a single positive integer.")
  }
  if (!rlang::is_logical(pad_na, n = 1)) {
    cli::cli_abort("{.arg pad_na} must be a single logical value.")
  }

  if (is.data.frame(x) || is.matrix(x)) {
    n_rows <- nrow(x)
    if (n_rows <= 2 * trim_length) {
      cli::cli_abort(
        "{.arg trim_length} is too large; it would affect all rows from the data."
      )
    }

    if (pad_na) {
      x[1:trim_length, ] <- NA
      x[(n_rows - trim_length + 1):n_rows, ] <- NA
      return(x)
    } else {
      return(x[(trim_length + 1):(n_rows - trim_length), , drop = FALSE])
    }
  } else if (is.atomic(x) && is.vector(x)) {
    n_len <- length(x)
    if (n_len <= 2 * trim_length) {
      cli::cli_abort(
        "{.arg trim_length} is too large; it would affect all elements from the vector."
      )
    }

    if (pad_na) {
      x[1:trim_length] <- NA
      x[(n_len - trim_length + 1):n_len] <- NA
      return(x)
    } else {
      return(x[(trim_length + 1):(n_len - trim_length)])
    }
  } else {
    cli::cli_abort("Input {.arg x} must be a vector, matrix, or data frame.")
  }
}

#' Downsample a Time Series Signal via Rolling Aggregation
#'
#' Reduces the sampling rate of a continuous time series by applying a rolling
#' aggregation function (median or mean) across non-overlapping windows.
#'
#' @details
#' **When to use this function versus `aggregate_by_time()`:**
#'
#' * Use `downsample_signal()` when working with a single, continuous numeric
#'   vector that has guaranteed regular intervals and no missing frames. That
#'   function relies on matrix reshaping and vector math, making it exceptionally
#'   fast for clean, pre-processed data.
#' * Use `aggregate_by_time()` when working with raw behavioral tracking data
#'   (e.g., OpenFace output) that may contain irregular timestamps, dropped
#'   frames, or missing rows. By binning based on the actual time variable,
#'   this function preserves the true chronological structure of the data and
#'   correctly leaves gaps where tracking was lost. It is also ideal for
#'   processing multiple numeric columns simultaneously.
#'
#' @param x A numeric vector representing the continuous time series signal.
#' @param factor A single positive integer indicating the downsampling factor.
#'   For example, a factor of 6 reduces 30Hz data to 5Hz.
#' @param method A character string specifying the aggregation method:
#'   "median" or "mean". Default is "median", which is highly robust to
#'   single-frame tracking glitches.
#' @param na.rm A logical indicating whether to remove missing values during
#'   aggregation. Default is `TRUE`.
#' @return A numeric vector representing the downsampled time series.
#' @export
downsample_signal <- function(
  x,
  factor,
  method = c("median", "mean"),
  na.rm = TRUE
) {
  if (!is.numeric(x)) {
    cli::cli_abort("{.arg x} must be a numeric vector.")
  }
  if (!rlang::is_integerish(factor, n = 1) || factor <= 1) {
    cli::cli_abort("{.arg factor} must be a single integer greater than 1.")
  }
  if (!rlang::is_logical(na.rm, n = 1)) {
    cli::cli_abort("{.arg na.rm} must be a single logical value.")
  }

  method <- match.arg(method)

  # Determine the number of complete windows
  n_windows <- floor(length(x) / factor)

  if (n_windows == 0) {
    cli::cli_abort(
      "The length of {.arg x} is smaller than the downsampling {.arg factor}."
    )
  }

  # Truncate the vector to fit perfectly into the windows
  x_truncated <- x[1:(n_windows * factor)]

  # Reshape into a matrix where each column is a window
  x_matrix <- matrix(x_truncated, nrow = factor, ncol = n_windows)

  # Apply the aggregation
  if (method == "median") {
    res <- apply(x_matrix, 2, stats::median, na.rm = na.rm)
  } else if (method == "mean") {
    res <- colMeans(x_matrix, na.rm = na.rm)
  }

  return(as.numeric(res))
}

# Main Functions ----------------------------------------------------------

#' Windowed Cross-Correlation
#'
#' Conduct a windowed cross-correlation analysis
#'
#' @param x A numeric vector containing a time series (same length as `y`).
#' @param y A numeric vector containing a time series (same length as `x`).
#' @param time An optional numeric vector representing the timestamps for the
#'   data. Must be the same length as `x` and `y`. If provided, the rolling
#'   window indices will be mapped directly to these timestamps in the results,
#'   which is highly recommended to maintain accurate timelines if edge
#'   artifacts were trimmed prior to analysis. Default is `NULL`.
#' @param window_size A positive integer indicating the size of each window,
#'   i.e., the number of elements in each window vector. Boker et al. suggest
#'   setting the window small enough so that the assumption can be made of
#'   little change in lead-lag relationships within the number of samples in the
#'   window but not so small that the reliability for the correlation estimate
#'   for each sample will be reduced.
#' @param lag_max A positive integer indicating the maximum lag to try between
#'   `x` and `y` windows. Boker et al. recommend selecting the greatest interval
#'   of time separating a behavior from participant `x` and a behavior from
#'   participant `y` that would be considered to be of interest.
#' @param window_increment A positive integer indicating the number of samples
#'   between successive changes in the window for the `x` vector. Can be made
#'   larger than 1 to reduce the number of rows in the output matrix. Boker et
#'   al. recommend setting the window increment as long as possible, but not so
#'   long that the relation between successive rows in the results matrix is
#'   lost. (default = `1`)
#' @param lag_increment A positive integer indicating the number of samples
#'   between successive changes in the window for the `y` vector (and thus also
#'   the interval of time separating successive columns in the results matrix).
#'   Boker et al. recommend setting the lag increment to the longest lag
#'   increment that still results in related change between successive columns.
#'   (default = `1`)
#' @param na.rm A logical indicating whether to remove missing values from the
#'   windows when calculating windowed cross-correlations. (default = `TRUE`)
#' @return A list object of class "wcc" containing the results matrix and useful
#'   summaries of it.
#' @export
wcc <- function(
  x,
  y,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1,
  na.rm = TRUE
) {
  # Assertions
  if (!is.numeric(x)) {
    cli::cli_abort("{.arg x} must be a numeric vector.")
  }
  if (!is.numeric(y)) {
    cli::cli_abort("{.arg y} must be a numeric vector.")
  }
  if (length(x) != length(y)) {
    cli::cli_abort("{.arg x} and {.arg y} must be the same length.")
  }

  if (!is.null(time)) {
    if (!is.numeric(time)) {
      cli::cli_abort("{.arg time} must be a numeric vector.")
    }
    if (length(time) != length(x)) {
      cli::cli_abort("{.arg time} must be the same length as {.arg x}.")
    }
  }

  if (!rlang::is_integerish(window_size, n = 1) || window_size <= 0) {
    cli::cli_abort("{.arg window_size} must be a single positive integer.")
  }
  if (!rlang::is_integerish(lag_max, n = 1) || lag_max <= 0) {
    cli::cli_abort("{.arg lag_max} must be a single positive integer.")
  }
  if (!rlang::is_integerish(window_increment, n = 1) || window_increment <= 0) {
    cli::cli_abort("{.arg window_increment} must be a single positive integer.")
  }
  if (!rlang::is_integerish(lag_increment, n = 1) || lag_increment <= 0) {
    cli::cli_abort("{.arg lag_increment} must be a single positive integer.")
  }
  if (!rlang::is_logical(na.rm, n = 1)) {
    cli::cli_abort("{.arg na.rm} must be a single logical value.")
  }

  # Conversions
  x <- as.double(x)
  y <- as.double(y)

  settings <- list(
    window_size = window_size,
    window_increment = window_increment,
    lag_max = lag_max,
    lag_increment = lag_increment,
    na.rm = na.rm,
    has_time = !is.null(time)
  )

  results_df <- create_wcc_df(
    x = x,
    y = y,
    time = time,
    settings = settings
  )

  out <- list(
    results_df = results_df,
    fisher_z = fisher_z(results_df),
    settings = settings
  )

  wcc_res(out)
}

#' Suggest WCC Hyperparameters
#'
#' Calculates principled starting values for Windowed Cross-Correlation parameters
#' based on the sampling rate of the data and the theoretical timing of the behaviors.
#'
#' @param sample_rate A numeric value indicating the sampling rate in Hertz (frames per second).
#' @param event_duration_sec The expected duration of a single behavioral event in seconds.
#'   Default is 2 (typical for brief conversational gestures).
#' @param max_delay_sec The maximum plausible reaction time between participants in seconds.
#'   Default is 3.
#' @param overlap_pct The desired percentage of overlap between consecutive time windows.
#'   Default is 0.5 (50 percent overlap).
#' @return A list of recommended parameters ready to be passed to `wcc()`.
#' @export
suggest_wcc_params <- function(
  sample_rate,
  event_duration_sec = 2,
  max_delay_sec = 3,
  overlap_pct = 0.5
) {
  suggested_window <- round((event_duration_sec * 4) * sample_rate)
  suggested_lag <- round(max_delay_sec * sample_rate)

  if (suggested_lag > (suggested_window / 2)) {
    cli::cli_warn(c(
      "The requested {.arg max_delay_sec} is too large relative to the {.arg event_duration_sec}.",
      "i" = "Capping {.arg lag_max} at half the {.arg window_size} to preserve statistical reliability."
    ))
    suggested_lag <- floor(suggested_window / 2)
  }

  suggested_w_inc <- max(1, round(suggested_window * (1 - overlap_pct)))

  cli::cli_h1("Suggested WCC Parameters")
  cli::cli_dl(c(
    "window_size" = "{suggested_window} ({round(suggested_window / sample_rate, 1)} seconds)",
    "lag_max" = "{suggested_lag} ({round(suggested_lag / sample_rate, 1)} seconds)",
    "window_increment" = "{suggested_w_inc} ({overlap_pct * 100}% overlap)",
    "lag_increment" = "1"
  ))

  invisible(list(
    window_size = suggested_window,
    lag_max = suggested_lag,
    window_increment = suggested_w_inc,
    lag_increment = 1
  ))
}

# Constructors ------------------------------------------------------------

new_wcc_res <- function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class = c("wcc_res", class(x)))
}

wcc_res <- function(x = list()) {
  new_wcc_res(x)
}

# S3 Methods --------------------------------------------------------------

#' Print method for wcc_res objects
#'
#' @param x An object of class "wcc_res".
#' @param ... Additional arguments (not used).
#' @export
print.wcc_res <- function(x, ...) {
  s <- x$settings
  n_windows <- length(unique(x$results_df$i))
  n_lags <- length(unique(x$results_df$tau))

  cli::cli_h1("Windowed Cross-Correlation Analysis")

  cli::cli_dl(c(
    "Total Windows" = "{n_windows}",
    "Total Lags Tested" = "{n_lags}",
    "Window Size" = "{s$window_size}",
    "Max Lag" = "{s$lag_max}",
    "Overall Fisher's Z" = "{round(x$fisher_z, 4)}"
  ))

  invisible(x)
}

#' Summary method for wcc_res objects
#'
#' @param object An object of class "wcc_res".
#' @param ... Additional arguments (not used).
#' @export
summary.wcc_res <- function(object, ...) {
  print(object)

  cli::cli_h2("Cross-Correlation Value Distribution")
  wcc_vals <- object$results_df$wcc
  q_vals <- stats::quantile(
    wcc_vals,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    na.rm = TRUE
  )

  print(round(q_vals, 4))

  n_na <- sum(is.na(wcc_vals))
  if (n_na > 0) {
    cli::cli_alert_warning("{n_na} missing value{?s} (NA) detected.")
  }

  invisible(object)
}

# Internal Helpers --------------------------------------------------------

#' @noRd
create_wcc_df <- function(x, y, time = NULL, settings) {
  n_x <- length(x)
  w_max <- settings$window_size
  w_inc <- settings$window_increment
  tau_max <- settings$lag_max
  tau_inc <- settings$lag_increment

  lags <- seq(-tau_max, tau_max, by = tau_inc)

  n_r <- floor((n_x - w_max - 2 * tau_max) / w_inc)
  n_c <- length(lags)

  results_df <- base::expand.grid(row = 1:n_r, col = 1:n_c) |>
    dplyr::mutate(
      i = 1 + tau_max + (row - 1) * w_inc,
      tau = lags[col],
      wcc = calc_wcc_cpp(
        x = x,
        y = y,
        i_vals = i,
        tau_vals = tau,
        w_max = w_max
      )
    )

  if (!is.null(time)) {
    results_df$i <- time[results_df$i]
  }

  results_df
}

#' @noRd
r_to_z <- function(r) {
  r[r == -1] <- -0.99
  r[r == 1] <- 0.99
  z <- 0.5 * base::log((1 + r) / (1 - r))
  z
}

#' @noRd
fisher_z <- function(results_df) {
  rvec <- results_df$wcc
  zvec <- r_to_z(rvec)
  azvec <- base::abs(zvec)
  res <- base::mean(azvec, na.rm = TRUE)

  res
}

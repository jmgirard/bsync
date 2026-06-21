# Main Functions ----------------------------------------------------------

#' Windowed Dynamic Time Warping
#'
#' Conduct a windowed dynamic time warping (WDTW) analysis to find the optimal
#' alignment distance between sliding windows of two time series.
#'
#' @param x A numeric vector containing a time series (same length as `y`).
#' @param y A numeric vector containing a time series (same length as `x`).
#' @param time An optional numeric vector representing the timestamps for the
#'   data. Must be the same length as `x` and `y`. Default is `NULL`.
#' @param window_size A positive integer indicating the size of each window.
#' @param lag_max A positive integer indicating the maximum lag to try.
#' @param window_increment A positive integer indicating the window shift increment. (default = `1`)
#' @param lag_increment A positive integer indicating the lag shift increment. (default = `1`)
#' @param scale_method Character string specifying how to standardize the data.
#'   "global" standardizes the entire time series before analysis. "local"
#'   standardizes within each sliding window. "none" applies no scaling. (default = `"global"`)
#' @param distance_metric Character string specifying the local cost function.
#'   "L1" uses absolute difference (Manhattan). "L2" uses squared difference
#'   (Euclidean). (default = `"L2"`)
#' @return A list object of class "wdtw_res".
#' @export
wdtw <- function(
  x,
  y,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1,
  scale_method = c("global", "local", "none"),
  distance_metric = c("L2", "L1")
) {
  scale_method <- match.arg(scale_method)
  distance_metric <- match.arg(distance_metric)

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

  if (scale_method == "global") {
    x <- as.numeric(base::scale(x))
    y <- as.numeric(base::scale(y))
  }

  x <- as.double(x)
  y <- as.double(y)

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

  settings <- list(
    window_size = window_size,
    window_increment = window_increment,
    lag_max = lag_max,
    lag_increment = lag_increment,
    scale_method = scale_method,
    distance_metric = distance_metric,
    has_time = !is.null(time)
  )

  results_df <- create_wdtw_df(x = x, y = y, time = time, settings = settings)
  mean_dist <- base::mean(results_df$dtw_dist, na.rm = TRUE)

  out <- list(
    results_df = results_df,
    mean_distance = mean_dist,
    settings = settings
  )

  new_wdtw_res(out)
}

# Constructors ------------------------------------------------------------

new_wdtw_res <- function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class = c("wdtw_res", class(x)))
}

# S3 Methods --------------------------------------------------------------

#' Print method for wdtw_res objects
#'
#' @param x An object of class "wdtw_res".
#' @param ... Additional arguments (not used).
#' @export
print.wdtw_res <- function(x, ...) {
  s <- x$settings
  n_windows <- length(unique(x$results_df$i))
  n_lags <- length(unique(x$results_df$tau))

  cli::cli_h1("Windowed Dynamic Time Warping Analysis")

  cli::cli_dl(c(
    "Total Windows" = "{n_windows}",
    "Total Lags Tested" = "{n_lags}",
    "Window Size" = "{s$window_size}",
    "Max Lag" = "{s$lag_max}",
    "Scale Method" = "{s$scale_method}",
    "Distance Metric" = "{s$distance_metric}",
    "Overall Mean Distance" = "{round(x$mean_distance, 4)}"
  ))

  invisible(x)
}

# Internal Helpers --------------------------------------------------------

#' @noRd
create_wdtw_df <- function(x, y, time = NULL, settings) {
  n_x <- length(x)
  w_max <- settings$window_size
  w_inc <- settings$window_increment
  tau_max <- settings$lag_max
  tau_inc <- settings$lag_increment

  use_l2 <- settings$distance_metric == "L2"
  local_scale <- settings$scale_method == "local"

  lags <- seq(-tau_max, tau_max, by = tau_inc)

  n_r <- floor((n_x - w_max - 2 * tau_max) / w_inc)
  n_c <- length(lags)

  results_df <- base::expand.grid(row = 1:n_r, col = 1:n_c) |>
    dplyr::mutate(
      i = 1 + tau_max + (row - 1) * w_inc,
      tau = lags[col],
      dtw_dist = calc_wdtw_cpp(
        x = x,
        y = y,
        i_vals = i,
        tau_vals = tau,
        w_max = w_max,
        use_l2 = use_l2,
        local_scale = local_scale
      )
    )

  if (!is.null(time)) {
    results_df$i <- time[results_df$i]
  }

  results_df
}

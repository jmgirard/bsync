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
#' @examples
#' \donttest{
#' # Windowed dynamic time warping. DTW is O(window^2) per cell, so this
#' # example runs on a short subset; use the full series in real analyses.
#' wdtw_res <- wdtw(
#'   x = sim_dyad$x_A[1:600],
#'   y = sim_dyad$x_B[1:600],
#'   window_size = 96,
#'   lag_max = 10
#' )
#' wdtw_res
#' }
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

  validate_series(x, y, time)
  validate_window_params(
    window_size, window_increment,
    lag_max = lag_max, lag_increment = lag_increment
  )

  if (scale_method == "global") {
    x <- as.numeric(base::scale(x))
    y <- as.numeric(base::scale(y))
  }

  x <- as.double(x)
  y <- as.double(y)

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
    aggregate  = c(mean_distance = mean_dist),
    settings   = settings
  )

  new_bsync_surface(out, "wdtw_res")
}

# S3 Methods --------------------------------------------------------------

#' Print method for wdtw_res objects
#'
#' @param x An object of class "wdtw_res".
#' @param ... Additional arguments (not used).
#' @return Returns `x` invisibly.
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
    "Overall Mean Distance" = "{round(x$aggregate[['mean_distance']], 4)}"
  ))

  invisible(x)
}

# Internal Helpers --------------------------------------------------------

#' @noRd
create_wdtw_df <- function(x, y, time = NULL, settings) {
  grid <- build_surface_grid(
    n_x = length(x),
    window_size = settings$window_size,
    window_increment = settings$window_increment,
    lag_max = settings$lag_max,
    lag_increment = settings$lag_increment,
    lagged = TRUE
  )

  use_l2 <- settings$distance_metric == "L2"
  local_scale <- settings$scale_method == "local"

  results_df <- data.frame(
    i = grid$i_vals,
    tau = grid$tau_vals,
    dtw_dist = calc_wdtw_cpp(
      x           = x,
      y           = y,
      i_vals      = grid$i_vals,
      tau_vals    = grid$tau_vals,
      w_max       = grid$w_max,
      use_l2      = use_l2,
      local_scale = local_scale
    )
  )

  if (!is.null(time)) {
    results_df$i <- time[results_df$i]
  }

  results_df
}

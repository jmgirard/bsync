# Main Functions ----------------------------------------------------------

#' Windowed Granger Causality
#'
#' Conduct a rolling windowed Granger Causality analysis to determine dynamic
#' leader-follower relationships between two continuous time series.
#'
#' @param x A numeric vector containing a time series (same length as `y`).
#' @param y A numeric vector containing a time series (same length as `x`).
#' @param time An optional numeric vector representing the timestamps for the data.
#' @param window_size A positive integer indicating the size of each window.
#' @param ar_order A positive integer specifying the Autoregressive (AR) order.
#'   This represents the maximum number of lags included in the prediction model.
#' @param window_increment A positive integer indicating the window shift increment. (default = `1`)
#' @return A list object of class "wgranger_res".
#' @examples
#' # Windowed Granger causality (no lag axis; directional F-statistics)
#' wgc_res <- wgranger(
#'   x = sim_dyad$x_A,
#'   y = sim_dyad$x_B,
#'   window_size = 96
#' )
#' wgc_res
#' @export
wgranger <- function(
  x,
  y,
  time = NULL,
  window_size,
  ar_order = 1,
  window_increment = 1
) {
  validate_series(x, y, time)
  validate_window_params(
    window_size, window_increment,
    ar_order = ar_order
  )

  x <- as.double(x)
  y <- as.double(y)

  settings <- list(
    window_size = window_size,
    ar_order = ar_order,
    window_increment = window_increment,
    has_time = !is.null(time)
  )

  results_df <- create_wgranger_df(x, y, time, settings)

  out <- list(
    results_df = results_df,
    aggregate = c(
      f_xy = base::mean(results_df$f_xy, na.rm = TRUE),
      f_yx = base::mean(results_df$f_yx, na.rm = TRUE)
    ),
    settings = settings
  )

  new_bsync_surface(out, "wgranger_res")
}

# S3 Methods --------------------------------------------------------------

#' Print method for wgranger_res objects
#'
#' @param x An object of class "wgranger_res".
#' @param ... Additional arguments (not used).
#' @return Returns `x` invisibly.
#' @export
print.wgranger_res <- function(x, ...) {
  s <- x$settings
  n_windows <- nrow(x$results_df)

  cli::cli_h1("Windowed Granger Causality Analysis")

  cli::cli_dl(c(
    "Total Windows" = "{n_windows}",
    "Window Size" = "{s$window_size}",
    "AR Order (Lags)" = "{s$ar_order}"
  ))

  invisible(x)
}

#' Summary method for wgranger_res objects
#'
#' @param object An object of class "wgranger_res".
#' @param ... Additional arguments (not used).
#' @return Returns `object` invisibly.
#' @export
summary.wgranger_res <- function(object, ...) {
  print(object)

  df <- object$results_df
  sig_x_to_y <- sum(df$p_xy < 0.05, na.rm = TRUE)
  sig_y_to_x <- sum(df$p_yx < 0.05, na.rm = TRUE)
  valid_n <- sum(!is.na(df$p_xy))

  cli::cli_h2("Significance Summary (p < 0.05)")

  if (valid_n > 0) {
    cli::cli_bullets(c(
      "*" = "'x' significantly predicts 'y' in {sig_x_to_y} windows ({round(sig_x_to_y / valid_n * 100, 1)}%)",
      "*" = "'y' significantly predicts 'x' in {sig_y_to_x} windows ({round(sig_y_to_x / valid_n * 100, 1)}%)"
    ))
  } else {
    cli::cli_alert_warning(
      "No valid metrics computed. Check window size and degrees of freedom."
    )
  }

  invisible(object)
}

# Internal Helpers --------------------------------------------------------

#' @noRd
create_wgranger_df <- function(x, y, time = NULL, settings) {
  grid <- build_surface_grid(
    n_x = length(x),
    window_size = settings$window_size,
    window_increment = settings$window_increment,
    lagged = FALSE
  )

  stats_df <- calc_wgranger_cpp(x, y, grid$i_vals, grid$w_max, settings$ar_order)

  results_df <- data.frame(i = grid$i_vals)
  results_df <- cbind(results_df, stats_df)

  if (!is.null(time)) {
    results_df$i <- time[results_df$i]
  }

  results_df
}

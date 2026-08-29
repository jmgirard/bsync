#' Plot wphase_res object
#'
#' @param x An object of class "wphase_res".
#' @param time_step A numeric value specifying the duration of each index.
#'   If not 1, axes will be converted from raw indices to time units. Default is 1.
#' @param color_low Character string specifying the color for a PLV of 0. Default is "#F7F7F7" (Off-white).
#' @param color_high Character string specifying the color for a PLV of 1. Default is "#2166AC" (Deep Blue).
#' @param show_zero_lag Logical indicating whether to draw a vertical line at lag = 0. Default is `TRUE`.
#' @param zero_line_color Character string specifying the color of the zero-lag line. Default is "black".
#' @param ... Additional arguments (not used).
#' @return A `ggplot2` plot object.
#' @examples
#' res <- wphase(sim_dyad$z_A, sim_dyad$z_B, window_size = 96, lag_max = 10)
#' plot(res)
#' @export
plot.wphase_res <- function(
  x,
  time_step = 1,
  color_low = "#F7F7F7",
  color_high = "#2166AC",
  show_zero_lag = TRUE,
  zero_line_color = "black",
  ...
) {
  build_surface_heatmap(
    df = x$results_df,
    fill_col = "plv",
    fill_scale = ggplot2::scale_fill_gradient(
      low      = color_low,
      high     = color_high,
      limits   = c(0, 1),
      na.value = "grey80",
      name     = "PLV"
    ),
    has_time = isTRUE(x$settings$has_time),
    time_step = time_step,
    show_zero_lag = show_zero_lag,
    zero_line_color = zero_line_color
  )
}

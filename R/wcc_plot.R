#' Plot wcc_res object
#'
#' @param x An object of class "wcc_res".
#' @param time_step A numeric value specifying the duration of each index.
#'   If not 1, axes will be converted from raw indices to time units. Default is 1.
#' @param color_low Character string specifying the color for a correlation of -1. Default is "#B2182B" (Deep Red).
#' @param color_mid Character string specifying the color for a correlation of 0. Default is "#F7F7F7" (Off-white).
#' @param color_high Character string specifying the color for a correlation of 1. Default is "#2166AC" (Deep Blue).
#' @param show_zero_lag Logical indicating whether to draw a vertical line at lag = 0. Default is `TRUE`.
#' @param zero_line_color Character string specifying the color of the zero-lag line. Default is "black".
#' @param ... Additional arguments (not used).
#' @return A `ggplot2` plot object.
#' @export
plot.wcc_res <- function(
  x,
  time_step = 1,
  color_low = "#B2182B",
  color_mid = "#F7F7F7",
  color_high = "#2166AC",
  show_zero_lag = TRUE,
  zero_line_color = "black",
  ...
) {
  build_surface_heatmap(
    df             = x$results_df,
    fill_col       = "wcc",
    fill_scale     = ggplot2::scale_fill_gradient2(
      low      = color_low,
      mid      = color_mid,
      high     = color_high,
      midpoint = 0,
      limits   = c(-1, 1),
      na.value = "grey80",
      name     = "WCC (r)"
    ),
    has_time       = isTRUE(x$settings$has_time),
    time_step      = time_step,
    show_zero_lag  = show_zero_lag,
    zero_line_color = zero_line_color
  )
}

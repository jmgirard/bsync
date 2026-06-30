#' Plot wdtw_res object
#'
#' @param x An object of class "wdtw_res".
#' @param time_step A numeric value specifying the duration of each index.
#'   If not 1, axes will be converted from raw indices to time units. Default is 1.
#' @param show_zero_lag Logical indicating whether to draw a vertical line at lag = 0. Default is `TRUE`.
#' @param zero_line_color Character string specifying the color of the zero-lag line. Default is "black".
#' @param ... Additional arguments (not used).
#' @return A `ggplot2` plot object.
#' @export
plot.wdtw_res <- function(
  x,
  time_step = 1,
  show_zero_lag = TRUE,
  zero_line_color = "black",
  ...
) {
  build_surface_heatmap(
    df = x$results_df,
    fill_col = "dtw_dist",
    fill_scale = ggplot2::scale_fill_gradientn(
      colors    = grDevices::hcl.colors(100, "viridis", rev = TRUE),
      na.value  = "grey80",
      name      = "DTW Distance"
    ),
    has_time = isTRUE(x$settings$has_time),
    time_step = time_step,
    show_zero_lag = show_zero_lag,
    zero_line_color = zero_line_color
  )
}

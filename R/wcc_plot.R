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
  df <- x$results_df
  has_time <- isTRUE(x$settings$has_time)

  # Scale lag to time if a step size is provided
  if (time_step != 1) {
    df$tau <- df$tau * time_step
    x_label <- "Lag (\u03c4) in Seconds"
  } else {
    x_label <- "Lag (\u03c4) Index"
  }

  # Only scale elapsed time if it was not natively provided
  if (!has_time && time_step != 1) {
    df$i <- df$i * time_step
    y_label <- "Elapsed Time (Seconds)"
  } else if (has_time) {
    y_label <- "Elapsed Time"
  } else {
    y_label <- "Elapsed Time Window Index"
  }

  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = tau, y = i, fill = wcc)) +
    ggplot2::geom_tile(na.rm = TRUE) +
    ggplot2::scale_fill_gradient2(
      low = color_low,
      mid = color_mid,
      high = color_high,
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "grey80",
      name = "WCC (r)"
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.5
      )
    ) +
    ggplot2::labs(
      x = x_label,
      y = y_label
    )

  # Add the optional zero-lag reference line
  if (show_zero_lag) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = 0,
        color = zero_line_color,
        linetype = "dashed",
        alpha = 0.5,
        linewidth = 0.5
      )
  }

  p
}

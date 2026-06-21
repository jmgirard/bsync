#' Plot wdtw_res object
#'
#' @param x An object of class "wdtw_res".
#' @param time_step A numeric value specifying the duration of each index.
#'   If not 1, axes will be converted from raw indices to time units. Default is 1.
#' @param show_zero_lag Logical indicating whether to draw a vertical line at lag = 0. Default is `TRUE`.
#' @param zero_line_color Character string specifying the color of the zero-lag line. Default is "black".
#' @param ... Additional arguments (not used).
#' @export
plot.wdtw_res <- function(
  x,
  time_step = 1,
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

  p <- ggplot2::ggplot(
    data = df,
    ggplot2::aes(x = tau, y = i, fill = dtw_dist)
  ) +
    ggplot2::geom_tile(na.rm = TRUE) +
    ggplot2::scale_fill_gradientn(
      colors = grDevices::hcl.colors(100, "viridis", rev = TRUE),
      na.value = "grey80",
      name = "DTW Distance"
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

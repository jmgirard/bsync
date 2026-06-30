#' Plot Surface with Optima Overlay
#'
#' @param surface_obj An object of class "wcc_res" or "wdtw_res".
#' @param optima_df A data frame of class "wcc_optima" or "wdtw_optima".
#' @param time_step A numeric value specifying the duration of each index. Default is 1.
#' @param line_color Character string specifying the color of the connecting line. Default is "black".
#' @param point_fill Character string specifying the inner fill of the optima points. Default is "black".
#' @param point_stroke Character string specifying the outer outline of the optima points. Default is "white".
#' @param show_zero_lag Logical indicating whether to draw a vertical line at lag = 0. Default is `TRUE`.
#' @param zero_line_color Character string specifying the color of the zero-lag line. Default is "black".
#' @param ... Additional arguments passed to the underlying plot method (e.g., custom colors for WCC).
#' @return A `ggplot2` plot object.
#' @examples
#' wcc_res <- wcc(sim_dyad$x_A, sim_dyad$x_B, window_size = 96, lag_max = 10)
#' optima <- pick_optima(wcc_res, L_size = 9)
#' plot_optima_overlay(wcc_res, optima)
#' @export
plot_optima_overlay <- function(
  surface_obj,
  optima_df,
  time_step = 1,
  line_color = "black",
  point_fill = "black",
  point_stroke = "white",
  show_zero_lag = TRUE,
  zero_line_color = "black",
  ...
) {
  p_base <- plot(
    surface_obj,
    time_step = time_step,
    show_zero_lag = show_zero_lag,
    zero_line_color = zero_line_color,
    ...
  )

  has_time <- isTRUE(surface_obj$settings$has_time)

  # Scale optimum coordinates to match the base plot
  if (time_step != 1) {
    optima_df$optimum_lag <- optima_df$optimum_lag * time_step
  }

  if (!has_time && time_step != 1) {
    optima_df$i <- optima_df$i * time_step
  }

  # Set a dynamic title based on the input object class
  is_wcc <- inherits(surface_obj, "wcc_res")
  plot_title <- ifelse(
    is_wcc,
    "Windowed Cross-Correlation with Optima Overlay",
    "Windowed Dynamic Time Warping with Optima Overlay"
  )

  p_final <- p_base +
    # The connecting line (drawn first so points sit on top)
    ggplot2::geom_path(
      data = optima_df,
      ggplot2::aes(x = optimum_lag, y = i),
      color = line_color,
      alpha = 0.6,
      na.rm = TRUE,
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    # The points with an outline for contrast
    ggplot2::geom_point(
      data = optima_df,
      ggplot2::aes(x = optimum_lag, y = i),
      shape = 21,
      fill = point_fill,
      color = point_stroke,
      size = 1.5,
      stroke = 0.4,
      na.rm = TRUE,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      title = plot_title
    )

  p_final
}

# Shared surface heatmap scaffold (M5) ---------------------------------------
#
# build_surface_heatmap() factors the time_step / axis-label / theme /
# zero-lag-reference-line scaffold shared by plot.wcc_res and plot.wdtw_res.
# The caller supplies only the fill scale (and optionally the fill column name).
# plot.wgranger_res is line-based, not a heatmap, and retains its own plot.

#' Shared scaffold for WCC / WDTW surface heatmap plots
#'
#' @param df `results_df` from a `wcc_res` or `wdtw_res` object.
#' @param fill_col Character; column name to map to fill (e.g. `"wcc"`,
#'   `"dtw_dist"`).
#' @param fill_scale A `ggplot2` scale object (e.g. `scale_fill_gradient2()`).
#' @param has_time Logical; whether `i` is a real timestamp (`settings$has_time`).
#' @param time_step Numeric; if not 1, scales raw indices to seconds.
#' @param show_zero_lag Logical; draw a dashed reference line at lag = 0.
#' @param zero_line_color Character; color for the zero-lag line.
#' @return A `ggplot2` object.
#' @noRd
build_surface_heatmap <- function(df, fill_col, fill_scale,
                                  has_time = FALSE, time_step = 1,
                                  show_zero_lag = TRUE,
                                  zero_line_color = "black") {
  if (time_step != 1) {
    df$tau <- df$tau * time_step
    x_label <- "Lag (τ) in Seconds"
  } else {
    x_label <- "Lag (τ) Index"
  }

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
    ggplot2::aes(x = tau, y = i, fill = .data[[fill_col]])
  ) +
    ggplot2::geom_tile(na.rm = TRUE) +
    fill_scale +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "black", fill = NA, linewidth = 0.5
      )
    ) +
    ggplot2::labs(x = x_label, y = y_label)

  if (show_zero_lag) {
    p <- p + ggplot2::geom_vline(
      xintercept = 0,
      color = zero_line_color,
      linetype = "dashed",
      alpha = 0.5,
      linewidth = 0.5
    )
  }

  p
}

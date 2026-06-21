#' Plot wgranger_res object
#'
#' @param x An object of class "wgranger_res".
#' @param time_step A numeric value specifying the duration of each index. Default is 1.
#' @param metric Character string specifying which metric to plot. Options are
#'   `"F"` for F-Statistics or `"p"` for negative log10 p-values. Default is `"F"`.
#' @param smooth Logical indicating whether to apply loess smoothing to the lines. Default is `FALSE`.
#' @param ... Additional arguments (not used).
#' @export
plot.wgranger_res <- function(
  x,
  time_step = 1,
  metric = c("F", "p"),
  smooth = FALSE,
  ...
) {
  metric <- match.arg(metric)
  df <- x$results_df
  has_time <- isTRUE(x$settings$has_time)

  if (!has_time && time_step != 1) {
    df$i <- df$i * time_step
    x_label <- "Elapsed Time (Seconds)"
  } else if (has_time) {
    x_label <- "Elapsed Time"
  } else {
    x_label <- "Elapsed Time Window Index"
  }

  # Reshape data manually to avoid heavy package dependencies
  df_xy <- data.frame(
    i = df$i,
    stat = if (metric == "F") df$f_xy else -log10(df$p_xy),
    direction = "x predicts y"
  )

  df_yx <- data.frame(
    i = df$i,
    stat = if (metric == "F") df$f_yx else -log10(df$p_yx),
    direction = "y predicts x"
  )

  plot_df <- rbind(df_xy, df_yx)

  y_label <- ifelse(metric == "F", "F-Statistic", "-log10(p-value)")
  plot_title <- "Rolling Windowed Granger Causality"

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = i, y = stat, color = direction)
  ) +
    ggplot2::scale_color_manual(
      values = c("x predicts y" = "#2166AC", "y predicts x" = "#B2182B")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = plot_title,
      x = x_label,
      y = y_label,
      color = "Information Flow"
    )

  if (smooth) {
    p <- p +
      ggplot2::geom_smooth(
        method = "loess",
        se = FALSE,
        linewidth = 1,
        na.rm = TRUE
      )
  } else {
    p <- p + ggplot2::geom_line(alpha = 0.8, linewidth = 0.8, na.rm = TRUE)
  }

  # If plotting p-values, add a horizontal reference line for p = 0.05
  if (metric == "p") {
    sig_level <- -log10(0.05)
    p <- p +
      ggplot2::geom_hline(
        yintercept = sig_level,
        color = "black",
        linetype = "dashed",
        alpha = 0.6
      ) +
      ggplot2::annotate(
        "text",
        x = min(plot_df$i, na.rm = TRUE),
        y = sig_level + 0.1,
        label = "p = 0.05",
        hjust = 0,
        vjust = 0,
        size = 3
      )
  }

  p
}

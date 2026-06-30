# Specification-curve plot for bsync_multiverse objects (M6) ------------------
#
# Implements a Simonsohn-style specification curve with:
#   Top panel : ES sorted by magnitude, colored by significance
#   Bottom panel: choice dashboard (which analytical choices each spec used)
#
# Uses only ggplot2 (already in Imports) + grid (base R recommended package).
# No additional CRAN dependencies needed.

#' Plot a synchrony multiverse specification curve
#'
#' Draws a Simonsohn-style specification curve for a `bsync_multiverse` object.
#' The top panel shows effect sizes sorted from smallest to largest, with
#' significant cells (p < .05) highlighted. The bottom panel is a choice
#' dashboard showing which analytic choices each specification used.
#'
#' @param x A `bsync_multiverse` object.
#' @param sig_color Color for significant cells (p < .05). Default: `"#2166AC"`.
#' @param insig_color Color for non-significant cells. Default: `"grey60"`.
#' @param active_color Fill for active choice tiles. Default: `"#2166AC"`.
#' @param point_size Size of ES points. Default: `1.5`.
#' @param top_frac Fraction of plot height allocated to the ES panel.
#'   Default: `0.55`.
#' @param ... Additional arguments (not used).
#' @return Returns `x` invisibly; draws to the active graphics device.
#' @seealso [synchrony_multiverse()], [tidy.bsync_multiverse()],
#'   [glance.bsync_multiverse()]
#' @export
plot.bsync_multiverse <- function(
  x,
  sig_color = "#2166AC",
  insig_color = "grey60",
  active_color = "#2166AC",
  point_size = 1.5,
  top_frac = 0.55,
  ...
) {
  gd <- x$grid[!is.na(x$grid$es), ]
  if (nrow(gd) == 0) {
    cli::cli_abort("No valid cells to plot (all ES values are NA).")
  }

  # Sort by ES
  gd <- gd[order(gd$es), ]
  gd$spec_rank <- seq_len(nrow(gd))
  gd$significant <- gd$p < 0.05

  # ------------------------------------------------------------------
  # Top panel: ES by spec rank
  # ------------------------------------------------------------------
  top_plot <- ggplot2::ggplot(
    gd,
    ggplot2::aes(x = spec_rank, y = es, color = significant)
  ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", color = "grey50",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = insig_color, "TRUE" = sig_color),
      labels = c("FALSE" = "No", "TRUE" = "Yes"),
      name   = "p < .05"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Effect Size (ES)",
      title = paste0(
        "Synchrony Multiverse  --  ",
        x$robustness$n_significant, " / ", x$robustness$n_cells,
        " cells significant  |  Median ES = ",
        round(x$robustness$median_es, 2)
      )
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "right"
    )

  # ------------------------------------------------------------------
  # Bottom panel: choice dashboard
  # ------------------------------------------------------------------
  # Build a long data frame: one row per (spec, dim, value) combination,
  # flagging which (dim, value) pairs are active for each spec.
  choice_dims <- c(
    "window_sec", "lag_sec", "increment_pct",
    "surrogate_method", "statistic"
  )
  # Keep only dims present in gd and with > 1 unique non-NA value
  choice_dims <- choice_dims[
    vapply(choice_dims, function(d) {
      d %in% names(gd) && length(unique(gd[[d]][!is.na(gd[[d]])])) >= 1
    }, logical(1))
  ]

  dim_labels <- c(
    window_sec = "Window\n(sec)",
    lag_sec = "Max Lag\n(sec)",
    increment_pct = "Increment\n(pct)",
    surrogate_method = "Surrogate\nMethod",
    statistic = "Statistic"
  )

  dash_rows <- list()
  y_pos <- 0L

  for (d in rev(choice_dims)) { # rev so window_sec ends up on top
    uvals <- sort(unique(gd[[d]][!is.na(gd[[d]])]))
    for (v in uvals) {
      y_pos <- y_pos + 1L
      dash_rows[[paste(d, v, sep = "_")]] <- data.frame(
        spec_rank = gd$spec_rank,
        y_val = y_pos,
        active = as.character(gd[[d]]) == as.character(v),
        row_label = paste0(dim_labels[d], ": ", v),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(dash_rows) == 0) {
    # No variation in choices; bottom panel is trivial
    bot_plot <- ggplot2::ggplot() +
      ggplot2::annotate("text",
        x = 0.5, y = 0.5,
        label = "No choice variation", size = 4
      ) +
      ggplot2::theme_void()
  } else {
    dash_data <- do.call(rbind, dash_rows)

    # Build ordered y-axis labels (same order as y_val)
    label_map <- unique(dash_data[, c("y_val", "row_label")])
    label_map <- label_map[order(label_map$y_val), ]
    y_breaks <- label_map$y_val
    y_labs <- label_map$row_label

    bot_plot <- ggplot2::ggplot(
      dash_data,
      ggplot2::aes(x = spec_rank, y = y_val, fill = active)
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_manual(
        values = c("FALSE" = "grey90", "TRUE" = active_color),
        guide  = "none"
      ) +
      ggplot2::scale_y_continuous(
        breaks = y_breaks,
        labels = y_labs,
        expand = ggplot2::expansion(add = 0.5)
      ) +
      ggplot2::labs(
        x = "Specification (sorted by ES)",
        y = NULL
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(color = "grey80", linewidth = 0.2),
        axis.text.y        = ggplot2::element_text(size = 7)
      )
  }

  # ------------------------------------------------------------------
  # Combine using grid viewports (base R grid package)
  # ------------------------------------------------------------------
  top_grob <- ggplot2::ggplotGrob(top_plot)
  bot_grob <- ggplot2::ggplotGrob(bot_plot)

  # Align widths so panels line up vertically
  max_w <- grid::unit.pmax(top_grob$widths, bot_grob$widths)
  top_grob$widths <- max_w
  bot_grob$widths <- max_w

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(
      nrow    = 2,
      ncol    = 1,
      heights = grid::unit(c(top_frac, 1 - top_frac), "npc")
    )
  ))

  grid::pushViewport(grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L))
  grid::grid.draw(top_grob)
  grid::popViewport()

  grid::pushViewport(grid::viewport(layout.pos.row = 2L, layout.pos.col = 1L))
  grid::grid.draw(bot_grob)
  grid::popViewport()

  grid::popViewport()

  invisible(x)
}

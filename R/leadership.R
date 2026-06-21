#' Calculate Leadership Asymmetry Index
#'
#' Computes a rolling index of leader-follower asymmetry from optimally picked lags.
#' A value of 1 indicates 'x' leads entirely, -1 indicates 'y' leads entirely,
#' and 0 indicates equal leading or simultaneous behavior.
#'
#' @param optima_obj An object of class "wcc_optima" or "wdtw_optima".
#' @param epoch_size A positive integer specifying the number of windows to group
#'   together to calculate the local asymmetry ratio.
#' @param min_valid A positive integer specifying the minimum number of valid (non-NA)
#'   optima required in an epoch to compute the index.
#' @return A data frame containing the rolling asymmetry index.
#' @export
leadership_asymmetry <- function(optima_obj, epoch_size = 10, min_valid = 3) {
  if (!inherits(optima_obj, c("wcc_optima", "wdtw_optima"))) {
    cli::cli_abort(
      "Input {.arg optima_obj} must be a {.cls wcc_optima} or {.cls wdtw_optima} object."
    )
  }

  if (!rlang::is_integerish(epoch_size, n = 1) || epoch_size <= 0) {
    cli::cli_abort("{.arg epoch_size} must be a single positive integer.")
  }

  df <- as.data.frame(optima_obj)
  n_rows <- nrow(df)

  # Initialize output vectors
  asymmetry_index <- rep(NA_real_, n_rows)
  x_leads_count <- rep(NA_integer_, n_rows)
  y_leads_count <- rep(NA_integer_, n_rows)

  half_epoch <- floor(epoch_size / 2)

  for (row in 1:n_rows) {
    start_idx <- max(1, row - half_epoch)
    end_idx <- min(n_rows, row + half_epoch)

    local_lags <- df$optimum_lag[start_idx:end_idx]
    local_lags <- local_lags[!is.na(local_lags)]

    if (length(local_lags) >= min_valid) {
      # Positive lag means x leads, negative means y leads
      x_leads <- sum(local_lags > 0)
      y_leads <- sum(local_lags < 0)

      x_leads_count[row] <- x_leads
      y_leads_count[row] <- y_leads

      # Calculate bounded index [-1, 1]
      total_directional <- x_leads + y_leads
      if (total_directional > 0) {
        asymmetry_index[row] <- (x_leads - y_leads) / total_directional
      } else {
        asymmetry_index[row] <- 0
      }
    }
  }

  out_df <- data.frame(
    i = df$i,
    x_leads_n = x_leads_count,
    y_leads_n = y_leads_count,
    asymmetry_index = asymmetry_index
  )
  out <- structure(out_df, class = c("bsync_lai", "data.frame"))
  attr(out, "has_time") <- isTRUE(attr(optima_obj, "has_time"))

  return(out)
}

# S3 Methods for bsync_lai ------------------------------------------------

#' Print method for bsync_lai objects
#'
#' @param x An object of class "bsync_lai".
#' @param n An integer specifying how many rows to print. Default is 5.
#' @param ... Additional arguments (not used).
#' @export
print.bsync_lai <- function(x, n = 5, ...) {
  df <- as.data.frame(x)
  valid_idx <- !is.na(df$asymmetry_index)
  valid_n <- sum(valid_idx)

  cli::cli_h1("Leadership Asymmetry Index (LAI)")

  if (valid_n > 0) {
    mean_lai <- mean(df$asymmetry_index, na.rm = TRUE)

    # Determine overall leader based on mean LAI
    if (mean_lai > 0.05) {
      leader_str <- "x predominantly leads"
    } else if (mean_lai < -0.05) {
      leader_str <- "y predominantly leads"
    } else {
      leader_str <- "Balanced / No clear dominant leader"
    }

    cli::cli_dl(c(
      "Valid Epochs Computed" = "{valid_n} out of {nrow(df)}",
      "Overall Mean LAI" = "{round(mean_lai, 3)}",
      "Overall Dynamics" = "{leader_str}"
    ))

    cli::cli_text("Showing the first {min(n, nrow(df))} result{?s}:")
    print(utils::head(df, n), row.names = FALSE)
  } else {
    cli::cli_alert_warning(
      "No valid LAI values could be computed. Check epoch size and input optima."
    )
  }

  invisible(x)
}

#' Plot bsync_lai object
#'
#' @param x An object of class "bsync_lai".
#' @param time_step A numeric value specifying the duration of each index. Default is 1.
#' @param line_color Character string specifying the color of the LAI line. Default is "black".
#' @param smooth Logical indicating whether to add a loess smoothing line. Default is `FALSE`.
#' @param ... Additional arguments (not used).
#' @export
plot.bsync_lai <- function(
  x,
  time_step = 1,
  line_color = "#2166AC",
  smooth = FALSE,
  ...
) {
  df <- as.data.frame(x)
  has_time <- isTRUE(attr(x, "has_time"))

  # Only scale elapsed time if it was not natively provided
  if (!has_time && time_step != 1) {
    df$i <- df$i * time_step
    x_label <- "Elapsed Time (Seconds)"
  } else if (has_time) {
    x_label <- "Elapsed Time"
  } else {
    x_label <- "Elapsed Time Window Index"
  }

  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = i, y = asymmetry_index)) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      linetype = "dashed",
      alpha = 0.6
    ) +
    ggplot2::geom_step(
      color = line_color,
      alpha = 0.8,
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-1, 1),
      breaks = c(-1, -0.5, 0, 0.5, 1)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "Dynamic Leadership Asymmetry",
      subtitle = "Values > 0 indicate 'x' leads; Values < 0 indicate 'y' leads",
      x = x_label,
      y = "Asymmetry Index"
    )

  if (smooth) {
    p <- p +
      ggplot2::geom_smooth(
        method = "loess",
        se = FALSE,
        color = "#B2182B",
        linewidth = 1,
        na.rm = TRUE
      )
  }

  p
}

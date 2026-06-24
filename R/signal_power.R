#' Evaluate Signal Power and Suggest Downsampling Rate
#'
#' Calculates the Power Spectral Density (PSD) of continuous time series using
#' Welch's method. It determines the frequency below which a specified proportion
#' of the total signal power is captured and recommends an optimal integer
#' downsampling factor that yields a clean sampling rate.
#'
#' @param x A numeric vector, a list of numeric vectors, or a data frame. If a
#'   data frame is provided, all numeric columns will be evaluated.
#' @param sample_rate A single positive number indicating the sampling rate in Hertz.
#' @param threshold A single numeric value between 0 and 1 indicating the cumulative
#'   proportion of power to capture. Default is 0.95 (95 percent).
#' @param plot A logical indicating whether to return a cumulative power plot. Default is `TRUE`.
#' @return A list containing the calculated cutoff frequencies, the recommended
#'   integer downsampling factor, the resulting target frequency, and optionally
#'   a `ggplot` object.
#' @export
evaluate_signal_power <- function(
  x,
  sample_rate,
  threshold = 0.95,
  plot = TRUE
) {
  if (
    !is.numeric(sample_rate) || length(sample_rate) != 1 || sample_rate <= 0
  ) {
    cli::cli_abort("{.arg sample_rate} must be a single positive number.")
  }
  if (
    !is.numeric(threshold) ||
      length(threshold) != 1 ||
      threshold <= 0 ||
      threshold >= 1
  ) {
    cli::cli_abort(
      "{.arg threshold} must be a single numeric value between 0 and 1."
    )
  }
  if (!rlang::is_logical(plot, n = 1)) {
    cli::cli_abort("{.arg plot} must be a single logical value.")
  }

  # Standardize input into a named list of numeric vectors
  if (is.numeric(x)) {
    x_list <- list(Signal = x)
  } else if (is.data.frame(x)) {
    x_list <- as.list(dplyr::select(x, dplyr::where(is.numeric)))
    if (length(x_list) == 0) {
      cli::cli_abort("{.arg x} contains no numeric columns.")
    }
  } else if (is.list(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      cli::cli_abort(
        "All elements in the list {.arg x} must be numeric vectors."
      )
    }
    x_list <- x
    if (is.null(names(x_list))) {
      names(x_list) <- paste0("Signal_", seq_along(x_list))
    }
  } else {
    cli::cli_abort(
      "{.arg x} must be a numeric vector, a list of numeric vectors, or a data frame."
    )
  }

  # Internal helper to calculate PSD for a single vector
  calc_psd <- function(vec) {
    vec_clean <- as.numeric(stats::na.omit(vec))

    if (length(vec_clean) < sample_rate) {
      return(NULL)
    }

    psd_res <- gsignal::pwelch(vec_clean, fs = sample_rate)
    cum_power <- cumsum(psd_res$spec) / sum(psd_res$spec)
    cutoff_idx <- which(cum_power >= threshold)[1]

    list(
      freqs = psd_res$freq,
      cum_power = cum_power,
      cutoff = psd_res$freq[cutoff_idx]
    )
  }

  results <- lapply(x_list, calc_psd)
  results <- results[!vapply(results, is.null, logical(1))]

  if (length(results) == 0) {
    cli::cli_abort(
      "No signals contained enough non-missing data to calculate PSD."
    )
  }

  all_cutoffs <- vapply(results, function(res) res$cutoff, numeric(1))
  is_multi <- length(all_cutoffs) > 1

  if (is_multi) {
    final_cutoff <- stats::quantile(all_cutoffs, probs = 0.95)
    cli::cli_h1("Dataset-Level Signal Power Evaluation")
    cli::cli_text("Evaluated {length(all_cutoffs)} signals.")
    cli::cli_text("95th percentile of cutoffs is {round(final_cutoff, 2)} Hz.")
  } else {
    final_cutoff <- all_cutoffs[[1]]
    cli::cli_h1("Signal Power Evaluation")
    cli::cli_text(
      "{threshold * 100}% of signal power is captured below {round(final_cutoff, 2)} Hz."
    )
  }

  # Calculate safe downsampling targets
  theoretical_min <- final_cutoff * 2
  safe_target_rate <- theoretical_min * 1.10
  max_integer_factor <- floor(sample_rate / safe_target_rate)

  best_factor <- max_integer_factor

  # Search for a clean integer Hz outcome up to 3 steps back
  if (max_integer_factor > 1) {
    search_range <- max_integer_factor:max(1, max_integer_factor - 3)
    for (f in search_range) {
      if (sample_rate %% f == 0) {
        best_factor <- f
        break
      }
    }
  }

  if (best_factor < 1) {
    cli::cli_alert_warning(
      "The signal contains high frequencies requiring a sampling rate near or above the original {sample_rate} Hz. Downsampling is not recommended."
    )
    recommended_rate <- sample_rate
    downsample_factor <- 1
  } else {
    recommended_rate <- sample_rate / best_factor
    downsample_factor <- best_factor

    cli::cli_alert_success(
      "Theoretical minimum rate is {round(theoretical_min, 2)} Hz."
    )
    cli::cli_alert_success(
      "Recommended integer downsampling factor: {downsample_factor} (Resulting Target Rate: {recommended_rate} Hz)."
    )
  }

  out <- list(
    primary_cutoff_freq = unname(final_cutoff),
    theoretical_min_rate = unname(theoretical_min),
    recommended_downsample_factor = unname(downsample_factor),
    recommended_target_rate = unname(recommended_rate)
  )

  if (is_multi) {
    out$all_cutoffs <- all_cutoffs
    out$summary_stats <- summary(all_cutoffs)
  }

  if (plot) {
    plot_data <- do.call(
      rbind,
      lapply(names(results), function(nm) {
        data.frame(
          Signal = nm,
          Frequency = results[[nm]]$freqs,
          CumulativePower = results[[nm]]$cum_power
        )
      })
    )

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data$Frequency,
        y = .data$CumulativePower,
        group = .data$Signal
      )
    ) +
      ggplot2::geom_hline(
        yintercept = threshold,
        color = "gray50",
        linetype = "dashed"
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        limits = c(0, 1)
      ) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.5
        )
      ) +
      ggplot2::labs(
        title = ifelse(
          is_multi,
          "Cumulative Power of Multiple Time Series",
          "Cumulative Power of Time Series"
        ),
        x = "Frequency (Hz)",
        y = "Cumulative Proportion of Power"
      )

    if (is_multi) {
      p <- p +
        ggplot2::geom_line(color = "#2166AC", alpha = 0.15, linewidth = 0.5) +
        ggplot2::geom_vline(
          xintercept = final_cutoff,
          color = "#B2182B",
          linetype = "dashed",
          linewidth = 1
        ) +
        ggplot2::annotate(
          "text",
          x = final_cutoff + 0.2,
          y = 0.1,
          label = "95th Percentile Cutoff",
          color = "#B2182B",
          hjust = 0,
          size = 3.5
        )
    } else {
      p <- p +
        ggplot2::geom_line(color = "#2166AC", linewidth = 1) +
        ggplot2::geom_vline(
          xintercept = final_cutoff,
          color = "#B2182B",
          linetype = "dashed",
          linewidth = 1
        )
    }

    out$plot <- p
    print(p)
  }

  invisible(out)
}

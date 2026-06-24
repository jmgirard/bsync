#' Auto-Tune WCC Parameters for a Dataset
#'
#' Automatically determines the optimal Windowed Cross-Correlation parameters
#' for a multi-dyad dataset by combining Power Spectral Density (PSD) analysis
#' with a surrogate-driven grid search.
#'
#' @details
#' **Reproducibility and Parallelization:**
#' This function involves random sampling (selecting dyads and generating surrogate
#' data). For reproducible results, call `set.seed()` before running this function.
#' To speed up computation, ensure you have set a parallel backend using the `future`
#' package (e.g., `future::plan(future::multisession)`) prior to execution.
#'
#' @param dyad_list A list of data frames, where each data frame represents a
#'   dyad and contains two numeric columns (the two time series).
#' @param sample_rate A single positive number indicating the sampling rate in Hertz.
#' @param n_tune_dyads Integer. The number of dyads to sample for the tuning phase.
#'   Default is 10 to provide a robust sample without excessive computation time.
#'   If the dataset has fewer than this number, all dyads are used.
#' @param n_surrogates Integer. Number of surrogates to generate per test.
#'   Default is 100, which provides a stable enough standard deviation to calculate
#'   standardized effect sizes during tuning.
#' @param surrogate_method Character string. "phase" (default) uses phase randomization,
#'   which preserves the power spectrum and is ideal for continuous physiological data.
#'   "circular" shifts the time series, which is better for preserving local
#'   autocorrelation in behavioral data.
#' @param trim_odd Logical. If `TRUE` and `surrogate_method = "phase"`, automatically
#'   drops the final observation of any odd-length time series to allow the Fourier
#'   transform to execute. Default is `FALSE`.
#' @param increment_pct Numeric value between 0.01 and 1.0. Determines the step size
#'   between successive windows as a percentage of the window size. Default is 0.05.
#' @param window_multipliers A numeric vector. Multipliers applied to the baseline cycle
#'   length to generate the grid of window sizes. Default is `c(0.5, 1.0, 2.0)`.
#' @param lag_multipliers A numeric vector. Multipliers applied to the window size
#'   to generate the grid of maximum lags. Default is `c(0.5, 1.0, 2.0)`.
#' @param min_window_size Integer. The absolute minimum number of observations required
#'   in a window to calculate a stable correlation. Default is 60.
#' @param progress Logical. If `TRUE` (default), displays a dynamic progress bar
#'   in the console during the grid search.
#' @return A list containing the optimal parameters and the full tuning grid results.
#' @export
autotune_wcc <- function(
  dyad_list,
  sample_rate,
  n_tune_dyads = 10,
  n_surrogates = 100,
  surrogate_method = c("phase", "circular"),
  trim_odd = FALSE,
  increment_pct = 0.05,
  window_multipliers = c(0.5, 1.0, 2.0),
  lag_multipliers = c(0.5, 1.0, 2.0),
  min_window_size = 60,
  progress = TRUE
) {
  surrogate_method <- match.arg(surrogate_method)

  if (!rlang::is_logical(progress, n = 1)) {
    cli::cli_abort("{.arg progress} must be a single logical value.")
  }
  if (!rlang::is_logical(trim_odd, n = 1)) {
    cli::cli_abort("{.arg trim_odd} must be a single logical value.")
  }
  if (increment_pct <= 0 || increment_pct > 1) {
    cli::cli_abort(
      "{.arg increment_pct} must be strictly greater than 0 and less than or equal to 1."
    )
  }

  cli::cli_h1("Starting WCC Auto-Tuning")

  # 1. Safely sample dyads to prevent memory bloat and long compute times
  n_total <- length(dyad_list)
  if (n_tune_dyads >= n_total) {
    tune_sample <- dyad_list
    cli::cli_alert_info("Using all {n_total} dyads for tuning.")
  } else {
    tune_sample <- sample(dyad_list, n_tune_dyads)
    cli::cli_alert_info(
      "Sampled {n_tune_dyads} dyads from {n_total} for tuning."
    )
  }

  if (surrogate_method == "phase" && trim_odd) {
    odd_count <- sum(vapply(
      tune_sample,
      function(df) nrow(df) %% 2 != 0,
      logical(1)
    ))
    if (odd_count > 0) {
      cli::cli_alert_warning(
        "Trimming 1 observation from {odd_count} odd-length dyads to enable phase randomization."
      )
    }
  }

  # 2. Analyze signal power only on the sampled data
  cli::cli_alert_info("Step 1: Analyzing signal power...")
  sample_signals <- unlist(
    lapply(tune_sample, function(df) list(df[[1]], df[[2]])),
    recursive = FALSE
  )

  psd_res <- evaluate_signal_power(
    x = sample_signals,
    sample_rate = sample_rate,
    plot = FALSE
  )
  baseline_cycle_sec <- 1 / psd_res$primary_cutoff_freq
  baseline_window <- round(baseline_cycle_sec * sample_rate)

  # 3. Generate the search grid dynamically
  cli::cli_alert_info("Step 2: Generating parameter grid...")
  test_windows <- round(baseline_window * window_multipliers)
  test_windows <- test_windows[test_windows >= min_window_size]

  if (length(test_windows) == 0) {
    cli::cli_abort(
      "Calculated window sizes are too small for reliable estimates. Minimum {min_window_size} samples required."
    )
  }

  grid <- base::expand.grid(
    window_size = unique(test_windows),
    lag_multiplier = unique(lag_multipliers)
  )
  grid$lag_max <- round(grid$window_size * grid$lag_multiplier)
  grid$lag_max[grid$lag_max < 1] <- 1
  grid$mean_effect_size <- NA

  cli::cli_alert_info(
    "Step 3: Evaluating {nrow(grid)} parameter combinations via surrogates..."
  )

  if (progress) {
    total_iterations <- nrow(grid) * length(tune_sample)
    cli::cli_progress_bar("Running Grid Search", total = total_iterations)
  }

  # 4. Grid Search
  for (i in seq_len(nrow(grid))) {
    w_size <- grid$window_size[i]
    l_max <- grid$lag_max[i]
    w_inc <- max(1, round(w_size * increment_pct))

    dyad_effects <- numeric(length(tune_sample))

    for (d in seq_along(tune_sample)) {
      x <- tune_sample[[d]][[1]]
      y <- tune_sample[[d]][[2]]

      if (surrogate_method == "phase" && trim_odd && length(y) %% 2 != 0) {
        x <- x[-length(x)]
        y <- y[-length(y)]
      }

      if (surrogate_method == "phase") {
        y_surrs <- generate_surrogate_phase(
          y,
          n_surrogates = n_surrogates,
          trim_odd = trim_odd
        )
      } else {
        y_surrs <- generate_surrogate_circular(
          y,
          n_surrogates = n_surrogates,
          lag_max = l_max
        )
      }

      surr_res <- wcc_surrogate(
        x = x,
        y = y,
        y_surrogates = y_surrs,
        window_size = w_size,
        lag_max = l_max,
        window_increment = w_inc
      )

      null_mean <- base::mean(surr_res$surrogate_z, na.rm = TRUE)
      null_sd <- stats::sd(surr_res$surrogate_z, na.rm = TRUE)

      if (is.na(null_sd) || null_sd == 0) {
        dyad_effects[d] <- 0
      } else {
        dyad_effects[d] <- (surr_res$observed_z - null_mean) / null_sd
      }

      if (progress) {
        cli::cli_progress_update()
      }
    }

    grid$mean_effect_size[i] <- base::mean(dyad_effects, na.rm = TRUE)
  }

  if (progress) {
    cli::cli_progress_done()
  }

  # 5. Extract Optimal Parameters
  optimal_params <- grid[which.max(grid$mean_effect_size), ]

  cli::cli_alert_success("Optimization complete.")
  cli::cli_dl(c(
    "Optimal Window Size" = "{optimal_params$window_size} ({optimal_params$window_size / sample_rate} sec)",
    "Optimal Max Lag" = "{optimal_params$lag_max} ({optimal_params$lag_max / sample_rate} sec)",
    "Maximized Effect Size (Z)" = "{round(optimal_params$mean_effect_size, 4)}"
  ))

  invisible(list(
    recommended_window_size = optimal_params$window_size,
    recommended_lag_max = optimal_params$lag_max,
    tuning_grid_results = grid
  ))
}

# Main Functions ----------------------------------------------------------

#' Windowed Cross-Correlation
#'
#' Conduct a windowed cross-correlation analysis
#'
#' @param x A numeric vector containing a time series (same length as `y`).
#' @param y A numeric vector containing a time series (same length as `x`).
#' @param time An optional numeric vector representing the timestamps for the
#'   data. Must be the same length as `x` and `y`. If provided, the rolling
#'   window indices will be mapped directly to these timestamps in the results,
#'   which is highly recommended to maintain accurate timelines if edge
#'   artifacts were trimmed prior to analysis. Default is `NULL`.
#' @param window_size A positive integer indicating the size of each window,
#'   i.e., the number of elements in each window vector. Boker et al. suggest
#'   setting the window small enough so that the assumption can be made of
#'   little change in lead-lag relationships within the number of samples in the
#'   window but not so small that the reliability for the correlation estimate
#'   for each sample will be reduced.
#' @param lag_max A positive integer indicating the maximum lag to try between
#'   `x` and `y` windows. Boker et al. recommend selecting the greatest interval
#'   of time separating a behavior from participant `x` and a behavior from
#'   participant `y` that would be considered to be of interest.
#' @param window_increment A positive integer indicating the number of samples
#'   between successive changes in the window for the `x` vector. Can be made
#'   larger than 1 to reduce the number of rows in the output matrix. Boker et
#'   al. recommend setting the window increment as long as possible, but not so
#'   long that the relation between successive rows in the results matrix is
#'   lost. (default = `1`)
#' @param lag_increment A positive integer indicating the number of samples
#'   between successive changes in the window for the `y` vector (and thus also
#'   the interval of time separating successive columns in the results matrix).
#'   Boker et al. recommend setting the lag increment to the longest lag
#'   increment that still results in related change between successive columns.
#'   (default = `1`)
#' @param na.rm A logical indicating whether to remove missing values from the
#'   windows when calculating windowed cross-correlations. (default = `TRUE`)
#' @param statistic A character string specifying how to aggregate the WCC
#'   surface into a single number. `"mean_abs_z"` (default) takes the mean of
#'   absolute Fisher's Z values over **all** windows and lags -- the SUSY *mean
#'   absolute Z* (Tschacher & Meier, 2020). `"peak"` takes the maximum absolute
#'   Fisher's Z across lags **within each window**, then averages those per-window
#'   peaks -- the rMEA *best-lag* convention (Boker et al., 2002). Both are
#'   larger-is-more-synchrony quantities. Pass the same value to
#'   `wcc_surrogate()` so the null distribution matches (see Invariant 2).
#' @return A list object of class "wcc_res" containing the results matrix and
#'   useful summaries of it.
#' @export
wcc <- function(
  x,
  y,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1,
  na.rm = TRUE,
  statistic = c("mean_abs_z", "peak")
) {
  statistic <- match.arg(statistic)

  validate_series(x, y, time)
  validate_window_params(
    window_size, window_increment,
    lag_max = lag_max, lag_increment = lag_increment,
    na.rm = na.rm
  )

  x <- as.double(x)
  y <- as.double(y)

  settings <- list(
    window_size = window_size,
    window_increment = window_increment,
    lag_max = lag_max,
    lag_increment = lag_increment,
    na.rm = na.rm,
    statistic = statistic,
    has_time = !is.null(time)
  )

  results_df <- create_wcc_df(
    x = x,
    y = y,
    time = time,
    settings = settings
  )

  agg_val <- wcc_aggregate(
    z = r_to_z(results_df$wcc),
    window_id = results_df$i,
    statistic = statistic
  )

  out <- list(
    results_df = results_df,
    aggregate  = stats::setNames(agg_val, statistic),
    settings   = settings
  )

  new_bsync_surface(out, "wcc_res")
}

#' Suggest WCC Hyperparameters
#'
#' Derives principled starting values for Windowed Cross-Correlation parameters
#' from the **measured signal** via PSD-based dominant timescale estimation.
#'
#' @details
#' **Dominant timescale.** When `event_duration_sec` is `NULL` (default), the
#' function estimates the dominant behavioral cycle from the measured signal via
#' [evaluate_signal_power()]: `event_duration_sec = 1 / primary_cutoff_freq`.
#' Pass a numeric value to override with your own theoretical estimate.
#'
#' **Window size (4-cycles heuristic).** `window_size = round(event_duration_sec
#' * 4 * sample_rate)` (Boker et al., 2002). Four cycles yields a stable
#' within-window correlation across the range of lead-lag relationships; two
#' cycles (the Nyquist minimum) is too noisy.
#'
#' **Hard constraints applied and reported:**
#' \itemize{
#'   \item `lag_max <= floor(window_size / 2)` -- the SUSY reliability constraint
#'     (`segment >= 2*maxlag`); beyond this the lagged windows share fewer than
#'     half their samples.
#'   \item `window_size <= floor(series_length / 2)` -- ensures at least two
#'     non-overlapping windows fit in the series.
#'   \item `window_size >= min_window_samples` -- a minimum-samples floor for a
#'     stable correlation estimate.
#' }
#' All violations produce informative messages, not silent changes.
#'
#' @param x Numeric vector; the reference time series.
#' @param y Numeric vector; the query time series (same length as `x`).
#' @param sample_rate A numeric value indicating the sampling rate in Hertz.
#' @param event_duration_sec Optional numeric override for the dominant
#'   behavioral cycle duration in seconds. Default `NULL` derives it from the
#'   PSD of `x` and `y` via [evaluate_signal_power()].
#' @param max_delay_sec The maximum plausible reaction time between participants
#'   in seconds, used to set the initial `lag_max`. Default is `3`.
#' @param overlap_pct The desired proportion of overlap between consecutive
#'   windows (0-1). Default is `0.5` (50\% overlap).
#' @param min_window_samples Minimum number of samples required in a window for
#'   a stable correlation. Default is `20`.
#' @return A named list with `window_size`, `lag_max`, `window_increment`, and
#'   `lag_increment`, ready to pass to [wcc()].
#' @export
suggest_wcc_params <- function(
  x,
  y,
  sample_rate,
  event_duration_sec = NULL,
  max_delay_sec = 3,
  overlap_pct = 0.5,
  min_window_samples = 20L
) {
  validate_series(x, y)
  if (!is.numeric(sample_rate) || length(sample_rate) != 1 || sample_rate <= 0) {
    cli::cli_abort("{.arg sample_rate} must be a single positive number.")
  }
  if (!is.null(event_duration_sec)) {
    if (!is.numeric(event_duration_sec) || length(event_duration_sec) != 1 ||
      event_duration_sec <= 0) {
      cli::cli_abort(
        "{.arg event_duration_sec} must be a single positive number or {.val NULL}."
      )
    }
  }
  if (!is.numeric(max_delay_sec) || length(max_delay_sec) != 1 || max_delay_sec <= 0) {
    cli::cli_abort("{.arg max_delay_sec} must be a single positive number.")
  }
  if (!is.numeric(overlap_pct) || length(overlap_pct) != 1 ||
    overlap_pct < 0 || overlap_pct >= 1) {
    cli::cli_abort("{.arg overlap_pct} must be a single number in [0, 1).")
  }
  if (!rlang::is_integerish(min_window_samples, n = 1) || min_window_samples < 2) {
    cli::cli_abort("{.arg min_window_samples} must be a single integer >= 2.")
  }

  n <- length(x)

  # --- Derive dominant timescale from PSD if not supplied ---------------
  if (is.null(event_duration_sec)) {
    psd_res <- evaluate_signal_power(
      x           = list(x = x, y = y),
      sample_rate = sample_rate,
      quiet       = TRUE
    )
    cycle_sec <- 1 / psd_res$primary_cutoff_freq
    cli::cli_alert_info(
      "PSD dominant cycle: {round(cycle_sec, 2)} s \\
       (cutoff {round(psd_res$primary_cutoff_freq, 2)} Hz)."
    )
  } else {
    cycle_sec <- event_duration_sec
    cli::cli_alert_info(
      "Using supplied {.arg event_duration_sec} = {event_duration_sec} s."
    )
  }

  # --- 4-cycles-per-window heuristic ------------------------------------
  suggested_window <- round(cycle_sec * 4 * sample_rate)

  # --- Enforce hard constraints with warnings ---------------------------
  # Constraint 1: minimum samples floor
  if (suggested_window < min_window_samples) {
    cli::cli_warn(c(
      "Derived window ({suggested_window} samples) is below {.arg min_window_samples} \\
       ({min_window_samples}).",
      "i" = "Increasing to {min_window_samples}."
    ))
    suggested_window <- as.integer(min_window_samples)
  }

  # Constraint 2: series-length ceiling (window <= series / 2)
  max_window <- floor(n / 2L)
  if (suggested_window > max_window) {
    cli::cli_warn(c(
      "Derived window ({suggested_window} samples) exceeds series_length/2 ({max_window}).",
      "i" = "Capping {.arg window_size} to {max_window}."
    ))
    suggested_window <- as.integer(max_window)
  }

  # Lag: convert seconds -> samples
  suggested_lag <- round(max_delay_sec * sample_rate)

  # Constraint 3: lag_max <= window / 2 (SUSY reliability constraint)
  max_lag <- floor(suggested_window / 2L)
  if (suggested_lag > max_lag) {
    cli::cli_warn(c(
      "Requested {.arg max_delay_sec} ({round(max_delay_sec, 2)} s = \\
       {round(suggested_lag)} samples) exceeds window_size/2 ({max_lag}).",
      "i" = "Capping {.arg lag_max} at {max_lag} (= {round(max_lag / sample_rate, 2)} s)."
    ))
    suggested_lag <- as.integer(max_lag)
  }

  suggested_w_inc <- max(1L, round(suggested_window * (1 - overlap_pct)))

  cli::cli_h1("Suggested WCC Parameters")
  cli::cli_dl(c(
    "window_size"      = "{suggested_window} ({round(suggested_window / sample_rate, 2)} s)",
    "lag_max"          = "{suggested_lag} ({round(suggested_lag / sample_rate, 2)} s)",
    "window_increment" = "{suggested_w_inc} ({overlap_pct * 100}% overlap)",
    "lag_increment"    = "1"
  ))

  invisible(list(
    window_size      = suggested_window,
    lag_max          = suggested_lag,
    window_increment = suggested_w_inc,
    lag_increment    = 1L
  ))
}

# S3 Methods --------------------------------------------------------------

#' Print method for wcc_res objects
#'
#' @param x An object of class "wcc_res".
#' @param ... Additional arguments (not used).
#' @return Returns `x` invisibly.
#' @export
print.wcc_res <- function(x, ...) {
  s <- x$settings
  n_windows <- length(unique(x$results_df$i))
  n_lags <- length(unique(x$results_df$tau))

  agg_label <- if (s$statistic == "peak") {
    "Mean Peak |Fisher's Z|"
  } else {
    "Mean |Fisher's Z|"
  }

  cli::cli_h1("Windowed Cross-Correlation Analysis")

  cli::cli_dl(c(
    "Total Windows" = "{n_windows}",
    "Total Lags Tested" = "{n_lags}",
    "Window Size" = "{s$window_size}",
    "Max Lag" = "{s$lag_max}",
    "{agg_label}" = "{round(x$aggregate[[1]], 4)}"
  ))

  invisible(x)
}

#' Summary method for wcc_res objects
#'
#' @param object An object of class "wcc_res".
#' @param ... Additional arguments (not used).
#' @return Returns `object` invisibly.
#' @export
summary.wcc_res <- function(object, ...) {
  print(object)

  cli::cli_h2("Cross-Correlation Value Distribution")
  wcc_vals <- object$results_df$wcc
  q_vals <- stats::quantile(
    wcc_vals,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    na.rm = TRUE
  )

  print(round(q_vals, 4))

  n_na <- sum(is.na(wcc_vals))
  if (n_na > 0) {
    cli::cli_alert_warning("{n_na} missing value{?s} (NA) detected.")
  }

  invisible(object)
}

# Internal Helpers --------------------------------------------------------

#' @noRd
create_wcc_df <- function(x, y, time = NULL, settings) {
  grid <- build_surface_grid(
    n_x = length(x),
    window_size = settings$window_size,
    window_increment = settings$window_increment,
    lag_max = settings$lag_max,
    lag_increment = settings$lag_increment,
    lagged = TRUE
  )

  results_df <- data.frame(
    i = grid$i_vals,
    tau = grid$tau_vals,
    wcc = calc_wcc_cpp(
      x        = x,
      y        = y,
      i_vals   = grid$i_vals,
      tau_vals = grid$tau_vals,
      w_max    = grid$w_max,
      na_rm    = settings$na.rm
    )
  )

  if (!is.null(time)) {
    results_df$i <- time[results_df$i]
  }

  results_df
}

#' @noRd
r_to_z <- function(r) {
  r_clamped <- pmax(pmin(r, 0.9999, na.rm = FALSE), -0.9999, na.rm = FALSE)
  z <- base::atanh(r_clamped)

  z
}

#' @noRd
wcc_aggregate <- function(z, window_id, statistic) {
  az <- base::abs(z)
  if (statistic == "mean_abs_z") {
    base::mean(az, na.rm = TRUE)
  } else {
    peaks <- tapply(az, window_id, function(v) {
      if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    })
    base::mean(peaks, na.rm = TRUE)
  }
}

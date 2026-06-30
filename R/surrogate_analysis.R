# =========================================================================
# === SURROGATE ANALYSIS WRAPPERS =========================================
# =========================================================================
# These functions evaluate the statistical significance of observed
# synchrony metrics against pre-generated matrices of surrogate data.
# =========================================================================

# -------------------------------------------------------------------------
# --- 1. Windowed Cross-Correlation (WCC) ---------------------------------
# -------------------------------------------------------------------------

#' Calculate Surrogate Windowed Cross-Correlations
#'
#' @details
#' The p-value is the proportion of surrogates whose aggregate statistic is **at least as
#' large as** the observed statistic. The aggregate — either `"mean_abs_z"` or `"peak"` —
#' is computed identically on the observed data and every surrogate via the same internal
#' helper, so the null distribution and the observed value are guaranteed to be directly
#' comparable (Invariant 2: surrogate nulls match the observed statistic).
#'
#' Pass the same `statistic` value you used in `wcc()` so that `observed_z` and the
#' surrogate draws use the same quantity.
#'
#' @param x A numeric vector containing a time series.
#' @param y A numeric vector containing a time series.
#' @param y_surrogates A matrix of surrogate time series for `y` (columns are surrogates).
#' @param time An optional numeric vector representing the timestamps for the data. Default is `NULL`.
#' @param window_size A positive integer indicating the size of each window.
#' @param lag_max A positive integer indicating the maximum lag to try.
#' @param window_increment A positive integer indicating the window shift increment. Default is 1.
#' @param lag_increment A positive integer indicating the lag shift increment. Default is 1.
#' @param na.rm A logical indicating whether to remove missing values. Default is `TRUE`.
#' @param statistic A character string specifying the aggregate statistic; must match the value
#'   passed to `wcc()`. `"mean_abs_z"` (default) or `"peak"`. See `wcc()` for details.
#' @return A list object of class "wcc_surr".
#' @export
wcc_surrogate <- function(
  x,
  y,
  y_surrogates,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1,
  na.rm = TRUE,
  statistic = c("mean_abs_z", "peak")
) {
  statistic <- match.arg(statistic)
  if (!is.matrix(y_surrogates)) {
    cli::cli_abort("{.arg y_surrogates} must be a matrix.")
  }
  if (nrow(y_surrogates) != length(y)) {
    cli::cli_abort(
      "{.arg y_surrogates} must have the same number of rows as length of {.arg y}."
    )
  }

  n_surrogates <- ncol(y_surrogates)

  # 1. Calculate observed WCC
  obs_wcc <- wcc(
    x = x,
    y = y,
    time = time,
    window_size = window_size,
    lag_max = lag_max,
    window_increment = window_increment,
    lag_increment = lag_increment,
    na.rm = na.rm,
    statistic = statistic
  )
  obs_z <- obs_wcc$aggregate[[1]]

  # 2. Build the structural grid ONCE for maximum speed
  grid <- build_surface_grid(
    n_x             = length(x),
    window_size      = window_size,
    window_increment = window_increment,
    lag_max          = lag_max,
    lag_increment    = lag_increment,
    lagged           = TRUE
  )
  i_vals   <- grid$i_vals
  tau_vals <- grid$tau_vals
  w_max_surr <- grid$w_max

  x_cpp <- as.double(x)

  # 3. Parallel-ready Loop using future.apply
  surrogate_zs <- future.apply::future_vapply(
    seq_len(n_surrogates),
    function(idx) {
      y_surr <- as.double(y_surrogates[, idx])

      wcc_vals <- calc_wcc_cpp(
        x = x_cpp,
        y = y_surr,
        i_vals = i_vals,
        tau_vals = tau_vals,
        w_max = w_max_surr,
        na_rm = na.rm
      )

      z_vals <- r_to_z(wcc_vals)

      # i_vals groups by window position (each unique i = one window)
      wcc_aggregate(z = z_vals, window_id = i_vals, statistic = statistic)
    },
    FUN.VALUE = numeric(1),
    future.seed = TRUE
  )

  p_val <- sum(surrogate_zs >= obs_z) / n_surrogates

  out <- list(
    observed_z = obs_z,
    surrogate_z = surrogate_zs,
    p_value = p_val,
    n_surrogates = n_surrogates,
    settings = obs_wcc$settings
  )

  structure(out, class = c("wcc_surr", "list"))
}

# -------------------------------------------------------------------------
# --- 2. Windowed Dynamic Time Warping (WDTW) -----------------------------
# -------------------------------------------------------------------------

#' Calculate Surrogate Windowed Dynamic Time Warping (WDTW)
#'
#' @details
#' The p-value is the proportion of surrogates whose aggregate statistic is **at most as
#' large as** the observed statistic (lower DTW distance = better alignment). The aggregate
#' is `mean(dtw_dist)` over all window × lag combinations — the same quantity stored in
#' `wdtw_res$mean_distance` — computed identically on both the observed data and every
#' surrogate, so the null distribution and the observed value are directly comparable.
#'
#' **`fast_method` warning:** when `fast_method = TRUE`, surrogates are evaluated at lag 0
#' only, while the observed statistic is computed over all lags. The null and observed
#' aggregates therefore cover different lag ranges, making the resulting p-value
#' approximate. Use only for quick exploratory checks, never for reporting.
#'
#' @param x A numeric vector containing the reference time series.
#' @param y A numeric vector containing the query time series.
#' @param y_surrogates A matrix of surrogate time series for `y` (columns are surrogates).
#' @param time An optional numeric vector representing the timestamps for the data. Default is `NULL`.
#' @param window_size A positive integer indicating the size of the rolling window.
#' @param lag_max A positive integer indicating the maximum lag to try.
#' @param window_increment A positive integer indicating the step size for the rolling window. Default is 1.
#' @param lag_increment A positive integer indicating the lag shift increment. Default is 1.
#' @param scale_method Character string specifying how to standardize the data. Default is `"global"`.
#' @param distance_metric Character string specifying the local cost function. Default is `"L2"`.
#' @param fast_method Logical. If `TRUE`, severely reduces computation time by only evaluating
#'   surrogate alignments at lag 0. **See Details for the statistical caveat.** Default is `FALSE`.
#' @return A list object of class "wdtw_surr".
#' @export
wdtw_surrogate <- function(
  x,
  y,
  y_surrogates,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1,
  scale_method = c("global", "local", "none"),
  distance_metric = c("L2", "L1"),
  fast_method = FALSE
) {
  scale_method <- match.arg(scale_method)
  distance_metric <- match.arg(distance_metric)
  use_l2 <- distance_metric == "L2"
  local_scale <- scale_method == "local"

  if (!is.matrix(y_surrogates)) {
    cli::cli_abort("{.arg y_surrogates} must be a matrix.")
  }
  if (nrow(y_surrogates) != length(y)) {
    cli::cli_abort(
      "{.arg y_surrogates} must have the same number of rows as length of {.arg y}."
    )
  }

  n_surrogates <- ncol(y_surrogates)

  # 1. Calculate observed WDTW
  obs_wdtw <- wdtw(
    x = x,
    y = y,
    time = time,
    window_size = window_size,
    lag_max = lag_max,
    window_increment = window_increment,
    lag_increment = lag_increment,
    scale_method = scale_method,
    distance_metric = distance_metric
  )
  obs_cost <- obs_wdtw$aggregate[["mean_distance"]]

  x_cpp <- as.double(x)
  if (scale_method == "global") {
    x_cpp <- as.numeric(base::scale(x))
  }

  # 2. Build grids for computation
  if (fast_method) {
    grid <- build_surface_grid(
      n_x             = length(x),
      window_size      = window_size,
      window_increment = window_increment,
      lag_max          = lag_max,
      lag_increment    = lag_increment,
      lagged           = FALSE
    )
    # Shift i_vals to account for lag_max offset (fast path uses lag=0 only)
    i_vals_surr   <- grid$i_vals + lag_max
    tau_vals_surr <- rep(0L, grid$n_r)
    cli::cli_alert_info(
      "Running fast method: Evaluating surrogates at lag 0 only."
    )
  } else {
    grid <- build_surface_grid(
      n_x             = length(x),
      window_size      = window_size,
      window_increment = window_increment,
      lag_max          = lag_max,
      lag_increment    = lag_increment,
      lagged           = TRUE
    )
    i_vals_surr   <- grid$i_vals
    tau_vals_surr <- grid$tau_vals
  }
  w_max_surr <- grid$w_max

  # 3. Parallel-ready loop using future.apply
  surrogate_costs <- future.apply::future_vapply(
    seq_len(n_surrogates),
    function(idx) {
      y_surr <- y_surrogates[, idx]

      if (scale_method == "global") {
        y_surr <- as.numeric(base::scale(y_surr))
      } else {
        y_surr <- as.double(y_surr)
      }

      surr_cost_vector <- calc_wdtw_cpp(
        x = x_cpp,
        y = y_surr,
        i_vals = i_vals_surr,
        tau_vals = tau_vals_surr,
        w_max = w_max_surr,
        use_l2 = use_l2,
        local_scale = local_scale
      )

      base::mean(surr_cost_vector, na.rm = TRUE)
    },
    FUN.VALUE = numeric(1),
    future.seed = TRUE
  )

  p_val <- sum(surrogate_costs <= obs_cost) / n_surrogates

  out <- list(
    observed_cost = obs_cost,
    surrogate_cost = surrogate_costs,
    p_value = p_val,
    n_surrogates = n_surrogates,
    settings = obs_wdtw$settings
  )

  structure(out, class = c("wdtw_surr", "list"))
}

# -------------------------------------------------------------------------
# --- 3. Windowed Granger Causality ---------------------------------------
# -------------------------------------------------------------------------

#' Calculate Surrogate Windowed Granger Causality
#'
#' @details
#' Two p-values are returned: one for x → y and one for y → x. Each is the proportion of
#' surrogates whose mean F-statistic across windows is **at least as large as** the
#' corresponding observed mean F-statistic. Both null distributions are built with the same
#' aggregate (`mean(f_xy)` and `mean(f_yx)`) as the observed statistics stored in
#' `wgranger_res$results_df`, so the null and observed values are directly comparable.
#'
#' @param x A numeric vector containing a time series.
#' @param y A numeric vector containing a time series.
#' @param y_surrogates A matrix of surrogate time series for `y` (columns are surrogates).
#' @param time An optional numeric vector representing the timestamps for the data. Default is `NULL`.
#' @param window_size A positive integer indicating the size of the rolling window.
#' @param ar_order A positive integer specifying the Autoregressive (AR) order. Default is 1.
#' @param window_increment A positive integer indicating the step size for the rolling window. Default is 1.
#' @return A list object of class "wgranger_surr".
#' @export
wgranger_surrogate <- function(
  x,
  y,
  y_surrogates,
  time = NULL,
  window_size,
  ar_order = 1,
  window_increment = 1
) {
  if (!is.matrix(y_surrogates)) {
    cli::cli_abort("{.arg y_surrogates} must be a matrix.")
  }
  if (nrow(y_surrogates) != length(y)) {
    cli::cli_abort(
      "{.arg y_surrogates} must have the same number of rows as length of {.arg y}."
    )
  }

  n_surrogates <- ncol(y_surrogates)

  # 1. Calculate observed Windowed Granger
  obs_wgranger <- wgranger(
    x = x,
    y = y,
    time = time,
    window_size = window_size,
    ar_order = ar_order,
    window_increment = window_increment
  )

  obs_f_xy <- obs_wgranger$aggregate[["f_xy"]]
  obs_f_yx <- obs_wgranger$aggregate[["f_yx"]]

  # 2. Setup structural grid
  grid <- build_surface_grid(
    n_x             = length(x),
    window_size      = window_size,
    window_increment = window_increment,
    lagged           = FALSE
  )
  i_vals     <- grid$i_vals
  w_max_surr <- grid$w_max

  x_cpp <- as.double(x)

  # 3. Parallel-ready loop using future.apply (returns 2 values per iteration)
  surr_matrix <- future.apply::future_vapply(
    seq_len(n_surrogates),
    function(idx) {
      y_surr <- as.double(y_surrogates[, idx])

      surr_stats <- calc_wgranger_cpp(
        x = x_cpp,
        y = y_surr,
        i_vals = i_vals,
        w_max = w_max_surr,
        p = ar_order
      )

      c(
        f_xy = base::mean(surr_stats$f_xy, na.rm = TRUE),
        f_yx = base::mean(surr_stats$f_yx, na.rm = TRUE)
      )
    },
    FUN.VALUE = numeric(2),
    future.seed = TRUE
  )

  # Extract the two rows into our separate numeric vectors
  surrogate_f_xy <- surr_matrix["f_xy", ]
  surrogate_f_yx <- surr_matrix["f_yx", ]

  # 4. Calculate empirical p-values (STANDARD LOGIC)
  p_val_xy <- sum(surrogate_f_xy >= obs_f_xy) / n_surrogates
  p_val_yx <- sum(surrogate_f_yx >= obs_f_yx) / n_surrogates

  out <- list(
    observed_f_xy = obs_f_xy,
    surrogate_f_xy = surrogate_f_xy,
    p_value_xy = p_val_xy,
    observed_f_yx = obs_f_yx,
    surrogate_f_yx = surrogate_f_yx,
    p_value_yx = p_val_yx,
    n_surrogates = n_surrogates,
    settings = obs_wgranger$settings
  )

  structure(out, class = c("wgranger_surr", "list"))
}

# =========================================================================
# === S3 PRINT METHODS ====================================================
# =========================================================================

#' Print method for wcc_surr objects
#'
#' @param x An object of class "wcc_surr".
#' @param ... Additional arguments (not used).
#' @return Returns `x` invisibly.
#' @export
print.wcc_surr <- function(x, ...) {
  cli::cli_h1("WCC Surrogate Analysis (Pseudo-Synchrony)")

  if (x$p_value == 0) {
    p_disp <- paste0("< ", 1 / x$n_surrogates)
  } else {
    p_disp <- as.character(round(x$p_value, 4))
  }

  agg_label <- if (x$settings$statistic == "peak") {
    "Mean Peak |Fisher's Z|"
  } else {
    "Mean |Fisher's Z|"
  }

  cli::cli_dl(c(
    "Permutations" = "{x$n_surrogates}",
    "Observed {agg_label}" = "{round(x$observed_z, 4)}",
    "Average Null {agg_label}" = "{round(mean(x$surrogate_z), 4)}",
    "Empirical p-value" = "{p_disp}"
  ))

  if (x$p_value < 0.05) {
    cli::cli_alert_success(
      "Observed synchrony is significantly greater than chance."
    )
  } else {
    cli::cli_alert_warning(
      "Observed synchrony is not significantly different from chance."
    )
  }

  if (x$n_surrogates < 1000) {
    cli::cli_alert_info(
      "Note: {x$n_surrogates} permutations may be too few for stable p-values.\n\tConsider setting `n_surrogates >= 1000` for final reporting."
    )
  }

  invisible(x)
}

#' Print method for wdtw_surr objects
#'
#' @param x An object of class "wdtw_surr".
#' @param ... Additional arguments (not used).
#' @return Returns `x` invisibly.
#' @export
print.wdtw_surr <- function(x, ...) {
  cli::cli_h1("WDTW Surrogate Analysis (Pseudo-Synchrony)")

  if (x$p_value == 0) {
    p_disp <- paste0("< ", 1 / x$n_surrogates)
  } else {
    p_disp <- as.character(round(x$p_value, 4))
  }

  cli::cli_dl(c(
    "Permutations" = "{x$n_surrogates}",
    "Observed Mean Cost" = "{round(x$observed_cost, 4)}",
    "Average Null Cost" = "{round(mean(x$surrogate_cost), 4)}",
    "Empirical p-value" = "{p_disp}"
  ))

  if (x$p_value < 0.05) {
    cli::cli_alert_success(
      "Observed cost is significantly lower than chance (stronger alignment)."
    )
  } else {
    cli::cli_alert_warning(
      "Observed cost is not significantly different from chance."
    )
  }

  if (x$n_surrogates < 1000) {
    cli::cli_alert_info(
      "Note: {x$n_surrogates} permutations may be too few for stable p-values.\n\tConsider setting `n_surrogates >= 1000` for final reporting."
    )
  }

  invisible(x)
}

#' Print method for wgranger_surr objects
#'
#' @param x An object of class "wgranger_surr".
#' @param ... Additional arguments (not used).
#' @return Returns `x` invisibly.
#' @export
print.wgranger_surr <- function(x, ...) {
  cli::cli_h1("Windowed Granger Surrogate Analysis")

  if (x$p_value_xy == 0) {
    p_disp_xy <- paste0("< ", 1 / x$n_surrogates)
  } else {
    p_disp_xy <- as.character(round(x$p_value_xy, 4))
  }
  if (x$p_value_yx == 0) {
    p_disp_yx <- paste0("< ", 1 / x$n_surrogates)
  } else {
    p_disp_yx <- as.character(round(x$p_value_yx, 4))
  }

  cli::cli_h2("Direction: x -> y")
  cli::cli_dl(c(
    "Permutations" = "{x$n_surrogates}",
    "Observed Mean F-statistic" = "{round(x$observed_f_xy, 4)}",
    "Average Null F-statistic" = "{round(mean(x$surrogate_f_xy), 4)}",
    "Empirical p-value" = "{p_disp_xy}"
  ))

  if (x$p_value_xy < 0.05) {
    cli::cli_alert_success(
      "Predictive power (x -> y) is significantly greater than chance."
    )
  } else {
    cli::cli_alert_warning(
      "Predictive power (x -> y) is not significantly different from chance."
    )
  }

  cli::cli_h2("Direction: y -> x")
  cli::cli_dl(c(
    "Observed Mean F-statistic" = "{round(x$observed_f_yx, 4)}",
    "Average Null F-statistic" = "{round(mean(x$surrogate_f_yx), 4)}",
    "Empirical p-value" = "{p_disp_yx}"
  ))

  if (x$p_value_yx < 0.05) {
    cli::cli_alert_success(
      "Predictive power (y -> x) is significantly greater than chance."
    )
  } else {
    cli::cli_alert_warning(
      "Predictive power (y -> x) is not significantly different from chance."
    )
  }

  if (x$n_surrogates < 1000) {
    cli::cli_alert_info(
      "Note: {x$n_surrogates} permutations may be too few for stable p-values.\n\tConsider setting `n_surrogates >= 1000` for final reporting."
    )
  }

  invisible(x)
}

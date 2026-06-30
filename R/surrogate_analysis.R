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

  # 2. Build the structural grid ONCE (hoistable in M6 multiverse)
  grid <- build_surface_grid(
    n_x = length(x),
    window_size = window_size,
    window_increment = window_increment,
    lag_max = lag_max,
    lag_increment = lag_increment,
    lagged = TRUE
  )

  x_cpp <- as.double(x)

  # aggregate-only compute_fn: core → aggregate, no results_df (Invariant 7)
  wcc_compute <- function(xv, y_col, g) {
    z_vals <- r_to_z(calc_wcc_cpp(
      x = xv, y = as.double(y_col),
      i_vals = g$i_vals, tau_vals = g$tau_vals,
      w_max = g$w_max, na_rm = na.rm
    ))
    wcc_aggregate(z = z_vals, window_id = g$i_vals, statistic = statistic)
  }

  # 3. Surrogate loop via shared engine
  surrogate_zs <- run_surrogate_engine(
    x = x_cpp, y_surrogates = y_surrogates, grid = grid,
    compute_fn = wcc_compute, fun_value = numeric(1)
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
#' `wdtw_res$aggregate[["mean_distance"]]` — computed identically on both the observed data and every
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

  # 2. Build grid (hoistable in M6 multiverse)
  grid <- build_surface_grid(
    n_x = length(x),
    window_size = window_size,
    window_increment = window_increment,
    lag_max = lag_max,
    lag_increment = lag_increment,
    lagged = TRUE
  )

  if (fast_method) {
    # Fast path: evaluate surrogates at lag 0 only, but over the *same* windows
    # the observed lagged surface uses (edges reserved at both ends). Take the
    # distinct window starts from the lagged grid — the col == 1 block, where
    # `row` varies fastest, so the first n_r entries are the unique positions —
    # and force tau = 0. (Using a lag-free grid here would shift windows past
    # the series end and over-count vs. the observed surface; see @details.)
    grid <- list(
      i_vals   = grid$i_vals[seq_len(grid$n_r)],
      tau_vals = rep(0L, grid$n_r),
      w_max    = grid$w_max,
      n_r      = grid$n_r
    )
    cli::cli_alert_info(
      "Running fast method: Evaluating surrogates at lag 0 only."
    )
  }

  # aggregate-only compute_fn (Invariant 7)
  wdtw_compute <- function(xv, y_col, g) {
    y_cpp_inner <- if (scale_method == "global") {
      as.numeric(base::scale(y_col))
    } else {
      as.double(y_col)
    }
    base::mean(calc_wdtw_cpp(
      x = xv, y = y_cpp_inner,
      i_vals = g$i_vals, tau_vals = g$tau_vals,
      w_max = g$w_max, use_l2 = use_l2, local_scale = local_scale
    ), na.rm = TRUE)
  }

  # 3. Surrogate loop via shared engine
  surrogate_costs <- run_surrogate_engine(
    x = x_cpp, y_surrogates = y_surrogates, grid = grid,
    compute_fn = wdtw_compute, fun_value = numeric(1)
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

  # 2. Setup structural grid (hoistable in M6 multiverse)
  grid <- build_surface_grid(
    n_x = length(x),
    window_size = window_size,
    window_increment = window_increment,
    lagged = FALSE
  )

  x_cpp <- as.double(x)

  # aggregate-only compute_fn returning named numeric(2) (Invariant 7)
  granger_compute <- function(xv, y_col, g) {
    surr_stats <- calc_wgranger_cpp(
      x = xv, y = as.double(y_col),
      i_vals = g$i_vals, w_max = g$w_max, p = ar_order
    )
    c(
      f_xy = base::mean(surr_stats$f_xy, na.rm = TRUE),
      f_yx = base::mean(surr_stats$f_yx, na.rm = TRUE)
    )
  }

  # 3. Surrogate loop via shared engine (returns 2 × n_surrogates matrix)
  surr_matrix <- run_surrogate_engine(
    x = x_cpp, y_surrogates = y_surrogates, grid = grid,
    compute_fn = granger_compute, fun_value = numeric(2)
  )

  surrogate_f_xy <- surr_matrix["f_xy", ]
  surrogate_f_yx <- surr_matrix["f_yx", ]

  # 4. Empirical p-values
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
    "Mean Peak Abs. Fisher's Z"
  } else {
    "Mean Abs. Fisher's Z"
  }

  cli::cli_dl(stats::setNames(
    c(
      "{x$n_surrogates}",
      "{round(x$observed_z, 4)}",
      "{round(mean(x$surrogate_z), 4)}",
      "{p_disp}"
    ),
    c(
      "Permutations",
      paste0("Observed ", agg_label),
      paste0("Average Null ", agg_label),
      "Empirical p-value"
    )
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

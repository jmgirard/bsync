# Synchrony Multiverse (M6) ---------------------------------------------------
#
# synchrony_multiverse() sweeps a seconds-specified parameter grid and evaluates
# synchrony under each specification using matched-null surrogates (Invariant 2).
#
# Key design points:
#  - Grid axes in seconds, converted per-cell to samples (window_size >= 1,
#    lag_max capped at floor(window_size/2))
#  - Surrogate reuse: one y_surrogates matrix per unique surrogate_method,
#    reused across every cell sharing it (efficiency seam from M5)
#  - ES polarity: WCC/Granger = upper-tail (higher better), WDTW = lower-tail
#  - Granger has two directional statistics (f_xy / f_yx); result columns
#    es_xy/p_xy/es_yx/p_yx are added when estimator = "wgranger"
#  - bsync_multiverse result is light (Invariant 7): tidy grid + settings +
#    robustness summary; no raw surrogate draws, no raw input stored

# select_specification() is the Phase B helper; defined in autotune.R


# synchrony_multiverse() -------------------------------------------------------

#' Synchrony Multiverse Analysis
#'
#' Sweeps a seconds-specified parameter grid across analytic choices and
#' evaluates each specification with a matched-null surrogate test (Invariant 2).
#' The headline metric is **effect size vs. the null** -- not raw synchrony, which
#' autocorrelation inflates.
#'
#' @details
#' **Grid construction.** Each vector argument (`window_sec`, `lag_sec`,
#' `increment_pct`, `statistic`, `surrogate_method`) is crossed into a full
#' parameter grid. Seconds are converted to samples per cell; `lag_max` is
#' hard-capped at `floor(window_size / 2)` to preserve statistical reliability.
#' Cells where the series is too short are silently skipped and appear as `NA`
#' rows in the output grid.
#'
#' **Surrogate reuse.** One surrogate matrix is generated per unique
#' `surrogate_method` and reused across every cell sharing that method; surrogate
#' cost does not multiply by grid size.
#'
#' **Effect size polarity.** For WCC and Granger, higher values indicate stronger
#' synchrony (upper-tail test): `ES = (obs - null_mean) / null_sd`. For WDTW,
#' lower distance is better (lower-tail test): `ES = (null_mean - obs) / null_sd`.
#' Both have the convention ES > 0 = evidence for synchrony.
#'
#' **Granger direction.** When `estimator = "wgranger"`, two sets of statistics
#' are returned: primary (`observed`, `null_mean`, `null_sd`, `es`, `p`) refer to
#' x -> y; additional columns `es_yx` and `p_yx` give y -> x.
#'
#' @param x Numeric vector; the reference time series.
#' @param y Numeric vector; the query time series (same length as `x`).
#' @param estimator Character string; windowed estimator: `"wcc"` (default),
#'   `"wdtw"`, or `"wgranger"`.
#' @param sample_rate Single positive number; sampling rate in Hz, used to
#'   convert `window_sec` and `lag_sec` to samples.
#' @param window_sec Numeric vector; window size(s) in seconds to sweep.
#' @param lag_sec Numeric vector; max lag(s) in seconds. Required for `"wcc"`
#'   and `"wdtw"`; ignored for `"wgranger"`. Capped at `window_sec / 2`
#'   per cell.
#' @param increment_pct Numeric vector; window increment(s) as a fraction of
#'   `window_size` (e.g., `0.1` = 10\% step). Default is `0.1`.
#' @param statistic Character vector; aggregate statistic for `"wcc"` only.
#'   One or both of `"mean_abs_z"` (SUSY) and `"peak"` (rMEA/Boker). Ignored
#'   for other estimators.
#' @param surrogate_method Character vector; surrogate generator(s): `"phase"`
#'   (preserves power spectrum) and/or `"circular"` (preserves autocorrelation).
#' @param n_surrogates Single positive integer; number of surrogates per cell.
#'   Default is `100`. Use >= 1000 for reporting.
#' @param ar_order Single positive integer; AR order for `"wgranger"`. Default
#'   is `1L`.
#' @param scale_method Character string; scaling for `"wdtw"`. Default is
#'   `"global"`.
#' @param distance_metric Character string; distance metric for `"wdtw"`.
#'   Default is `"L2"`.
#' @param na.rm Logical; passed to `"wcc"`. Default is `TRUE`.
#' @return A `bsync_multiverse` object with:
#'   \describe{
#'     \item{`$grid`}{[tibble::tibble()] with one row per parameter cell:
#'       specification columns, `window_size`/`lag_max`/`window_increment`
#'       (samples), `n_windows`, `observed`, `null_mean`, `null_sd`, `es`, `p`
#'       (plus `es_yx`/`p_yx` for Granger).}
#'     \item{`$settings`}{Named list of call-level inputs.}
#'     \item{`$robustness`}{Named list: `n_cells`, `n_significant`,
#'       `pct_significant`, `median_es`, `iqr_es`, `sign_consistent`
#'       (proportion of significant cells with ES > 0).}
#'   }
#' @seealso [autotune_wcc()], [suggest_wcc_params()], [plot.bsync_multiverse()],
#'   [tidy.bsync_multiverse()], [glance.bsync_multiverse()]
#' @export
synchrony_multiverse <- function(
  x, y,
  estimator = c("wcc", "wdtw", "wgranger"),
  sample_rate,
  window_sec,
  lag_sec = NULL,
  increment_pct = 0.1,
  statistic = "mean_abs_z",
  surrogate_method = "phase",
  n_surrogates = 100L,
  ar_order = 1L,
  scale_method = c("global", "local", "none"),
  distance_metric = c("L2", "L1"),
  na.rm = TRUE
) {
  estimator <- match.arg(estimator)
  scale_method <- match.arg(scale_method)
  distance_metric <- match.arg(distance_metric)
  statistic <- match.arg(statistic, choices = c("mean_abs_z", "peak"), several.ok = TRUE)
  surrogate_method <- match.arg(
    surrogate_method,
    choices = c("phase", "circular"), several.ok = TRUE
  )

  # --- Input validation -------------------------------------------------
  validate_series(x, y)
  if (!is.numeric(sample_rate) || length(sample_rate) != 1 || sample_rate <= 0) {
    cli::cli_abort("{.arg sample_rate} must be a single positive number.")
  }
  if (!is.numeric(window_sec) || length(window_sec) < 1 || any(window_sec <= 0)) {
    cli::cli_abort("{.arg window_sec} must be a positive numeric vector.")
  }
  if (!is.numeric(increment_pct) || any(increment_pct <= 0) || any(increment_pct > 1)) {
    cli::cli_abort("{.arg increment_pct} must be a numeric vector in (0, 1].")
  }
  if (!rlang::is_integerish(n_surrogates, n = 1) || n_surrogates < 1) {
    cli::cli_abort("{.arg n_surrogates} must be a single positive integer.")
  }

  is_granger <- estimator == "wgranger"

  if (!is_granger) {
    if (is.null(lag_sec)) {
      cli::cli_abort(
        "{.arg lag_sec} is required for {.arg estimator = {.val {estimator}}}."
      )
    }
    if (!is.numeric(lag_sec) || any(lag_sec <= 0)) {
      cli::cli_abort("{.arg lag_sec} must be a positive numeric vector.")
    }
  }

  n <- length(x)

  # --- Build parameter grid ---------------------------------------------
  if (is_granger) {
    grid_params <- base::expand.grid(
      window_sec       = window_sec,
      lag_sec          = NA_real_,
      increment_pct    = increment_pct,
      surrogate_method = surrogate_method,
      statistic        = NA_character_,
      stringsAsFactors = FALSE
    )
  } else if (estimator == "wcc") {
    grid_params <- base::expand.grid(
      window_sec       = window_sec,
      lag_sec          = lag_sec,
      increment_pct    = increment_pct,
      surrogate_method = surrogate_method,
      statistic        = statistic,
      stringsAsFactors = FALSE
    )
  } else {
    # wdtw
    grid_params <- base::expand.grid(
      window_sec       = window_sec,
      lag_sec          = lag_sec,
      increment_pct    = increment_pct,
      surrogate_method = surrogate_method,
      statistic        = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  # --- Pre-generate one surrogate matrix per surrogate_method -----------
  # For circular, we use the largest lag_max in the grid as the minimum shift
  # (more conservative; valid for all cells). For Granger with no lag_sec,
  # use n/4 as a safe default.
  unique_methods <- unique(grid_params$surrogate_method)

  if (!is_granger) {
    # Largest lag_max across the grid (in samples, pre-capped at window/2)
    max_lag_samp <- max(vapply(seq_len(nrow(grid_params)), function(ci) {
      w_samp <- max(1L, round(grid_params$window_sec[ci] * sample_rate))
      l_raw <- max(1L, round(grid_params$lag_sec[ci] * sample_rate))
      min(l_raw, floor(w_samp / 2L))
    }, numeric(1)))
    circ_lag_max <- max_lag_samp
  } else {
    circ_lag_max <- max(1L, round(n / 4L))
  }

  surr_matrices <- lapply(unique_methods, function(method) {
    if (method == "phase") {
      generate_surrogate_phase(y, n_surrogates = n_surrogates)
    } else {
      generate_surrogate_circular(y,
        n_surrogates = n_surrogates,
        lag_max = circ_lag_max
      )
    }
  })
  names(surr_matrices) <- unique_methods

  # --- Evaluate each cell -----------------------------------------------
  n_cells <- nrow(grid_params)

  window_size_samp <- integer(n_cells)
  lag_max_samp <- integer(n_cells)
  window_inc_samp <- integer(n_cells)
  n_windows_vec <- integer(n_cells)
  observed_vec <- numeric(n_cells)
  null_mean_vec <- numeric(n_cells)
  null_sd_vec <- numeric(n_cells)
  es_vec <- numeric(n_cells)
  p_vec <- numeric(n_cells)
  skipped_vec <- logical(n_cells)

  if (is_granger) {
    observed_yx_vec <- numeric(n_cells)
    null_mean_yx_vec <- numeric(n_cells)
    null_sd_yx_vec <- numeric(n_cells)
    es_yx_vec <- numeric(n_cells)
    p_yx_vec <- numeric(n_cells)
  }

  for (ci in seq_len(n_cells)) {
    w_sec <- grid_params$window_sec[ci]
    l_sec <- grid_params$lag_sec[ci]
    inc_pct <- grid_params$increment_pct[ci]
    s_method <- grid_params$surrogate_method[ci]
    stat <- grid_params$statistic[ci]

    # Convert seconds -> samples
    w_samp <- max(1L, round(w_sec * sample_rate))
    l_samp <- if (is_granger) {
      NA_integer_
    } else {
      l_raw <- max(1L, round(l_sec * sample_rate))
      as.integer(min(l_raw, floor(w_samp / 2L)))
    }
    inc_samp <- max(1L, round(w_samp * inc_pct))

    window_size_samp[ci] <- w_samp
    lag_max_samp[ci] <- if (is_granger) NA_integer_ else l_samp
    window_inc_samp[ci] <- inc_samp

    # Guard against short series (build_surface_grid would abort)
    min_n <- if (is_granger) {
      (w_samp - 1L) + inc_samp # at least one window
    } else {
      (w_samp - 1L) + 2L * l_samp + inc_samp
    }
    if (n < min_n) {
      skipped_vec[ci] <- TRUE
      observed_vec[ci] <- NA_real_
      null_mean_vec[ci] <- NA_real_
      null_sd_vec[ci] <- NA_real_
      es_vec[ci] <- NA_real_
      p_vec[ci] <- NA_real_
      n_windows_vec[ci] <- 0L
      if (is_granger) {
        observed_yx_vec[ci] <- NA_real_
        null_mean_yx_vec[ci] <- NA_real_
        null_sd_yx_vec[ci] <- NA_real_
        es_yx_vec[ci] <- NA_real_
        p_yx_vec[ci] <- NA_real_
      }
      next
    }

    y_surr <- surr_matrices[[s_method]]

    cell <- tryCatch(
      {
        if (estimator == "wcc") {
          .mv_wcc_cell(x, y, w_samp, l_samp, inc_samp, stat, na.rm, y_surr)
        } else if (estimator == "wdtw") {
          .mv_wdtw_cell(
            x, y, w_samp, l_samp, inc_samp, scale_method,
            distance_metric, y_surr
          )
        } else {
          .mv_granger_cell(x, y, w_samp, inc_samp, ar_order, y_surr)
        }
      },
      error = function(e) {
        list(
          observed = NA_real_, null_mean = NA_real_, null_sd = NA_real_,
          es = NA_real_, p = NA_real_, n_windows = 0L,
          observed_yx = NA_real_, null_mean_yx = NA_real_, null_sd_yx = NA_real_,
          es_yx = NA_real_, p_yx = NA_real_
        )
      }
    )

    n_windows_vec[ci] <- cell$n_windows
    observed_vec[ci] <- cell$observed
    null_mean_vec[ci] <- cell$null_mean
    null_sd_vec[ci] <- cell$null_sd
    es_vec[ci] <- cell$es
    p_vec[ci] <- cell$p

    if (is_granger) {
      observed_yx_vec[ci] <- cell$observed_yx
      null_mean_yx_vec[ci] <- cell$null_mean_yx
      null_sd_yx_vec[ci] <- cell$null_sd_yx
      es_yx_vec[ci] <- cell$es_yx
      p_yx_vec[ci] <- cell$p_yx
    }
  }

  # Warn about skipped cells
  n_skipped <- sum(skipped_vec)
  if (n_skipped > 0) {
    cli::cli_alert_warning(
      "{n_skipped} cell{?s} skipped: series too short for the requested parameters."
    )
  }

  # --- Assemble tidy grid -----------------------------------------------
  grid_tbl <- tibble::tibble(
    estimator        = estimator,
    window_sec       = grid_params$window_sec,
    lag_sec          = grid_params$lag_sec,
    increment_pct    = grid_params$increment_pct,
    surrogate_method = grid_params$surrogate_method,
    statistic        = grid_params$statistic,
    window_size      = window_size_samp,
    lag_max          = lag_max_samp,
    window_increment = window_inc_samp,
    n_windows        = n_windows_vec,
    observed         = observed_vec,
    null_mean        = null_mean_vec,
    null_sd          = null_sd_vec,
    es               = es_vec,
    p                = p_vec
  )

  if (is_granger) {
    grid_tbl$observed_yx <- observed_yx_vec
    grid_tbl$null_mean_yx <- null_mean_yx_vec
    grid_tbl$null_sd_yx <- null_sd_yx_vec
    grid_tbl$es_yx <- es_yx_vec
    grid_tbl$p_yx <- p_yx_vec
  }

  # --- Robustness summary -----------------------------------------------
  valid <- !skipped_vec & !is.na(es_vec)
  valid_es <- es_vec[valid]
  valid_p <- p_vec[valid]
  n_valid <- length(valid_es)
  n_sig <- sum(valid_p < 0.05, na.rm = TRUE)

  robustness <- list(
    n_cells = n_valid,
    n_significant = n_sig,
    pct_significant = if (n_valid > 0) n_sig / n_valid else NA_real_,
    median_es = if (n_valid > 0) stats::median(valid_es, na.rm = TRUE) else NA_real_,
    iqr_es = if (n_valid > 0) stats::IQR(valid_es, na.rm = TRUE) else NA_real_,
    sign_consistent = if (n_sig > 0) {
      mean(valid_es[valid_p < 0.05] > 0, na.rm = TRUE)
    } else {
      NA_real_
    }
  )

  settings <- list(
    estimator    = estimator,
    sample_rate  = sample_rate,
    n_surrogates = n_surrogates,
    n_cells      = n_cells
  )

  out <- list(grid = grid_tbl, settings = settings, robustness = robustness)
  structure(out, class = c("bsync_multiverse", "list"))
}


# Per-estimator cell helpers ---------------------------------------------------

#' @noRd
.mv_wcc_cell <- function(x, y, window_size, lag_max, window_increment, statistic,
                         na.rm, y_surr) {
  surr_res <- wcc_surrogate(
    x = x, y = y,
    y_surrogates = y_surr,
    window_size = window_size,
    lag_max = lag_max,
    window_increment = window_increment,
    na.rm = na.rm,
    statistic = statistic
  )
  obs <- surr_res$observed_z
  surr_vals <- surr_res$surrogate_z
  null_mean <- mean(surr_vals, na.rm = TRUE)
  null_sd <- stats::sd(surr_vals, na.rm = TRUE)
  es <- if (!is.na(null_sd) && null_sd > 0) (obs - null_mean) / null_sd else NA_real_
  p <- surr_res$p_value
  # n_windows via grid
  n_r <- build_surface_grid(
    n_x = length(x), window_size = window_size,
    window_increment = window_increment, lag_max = lag_max,
    lag_increment = 1L, lagged = TRUE
  )$n_r
  list(
    observed = obs, null_mean = null_mean, null_sd = null_sd,
    es = es, p = p, n_windows = n_r
  )
}

#' @noRd
.mv_wdtw_cell <- function(x, y, window_size, lag_max, window_increment,
                          scale_method, distance_metric, y_surr) {
  surr_res <- wdtw_surrogate(
    x = x, y = y,
    y_surrogates = y_surr,
    window_size = window_size,
    lag_max = lag_max,
    window_increment = window_increment,
    scale_method = scale_method,
    distance_metric = distance_metric
  )
  obs <- surr_res$observed_cost
  surr_vals <- surr_res$surrogate_cost
  null_mean <- mean(surr_vals, na.rm = TRUE)
  null_sd <- stats::sd(surr_vals, na.rm = TRUE)
  # WDTW: lower distance = better -> ES = (null_mean - obs) / null_sd
  es <- if (!is.na(null_sd) && null_sd > 0) (null_mean - obs) / null_sd else NA_real_
  # p = proportion of surrogates <= obs (lower tail)
  p <- surr_res$p_value
  n_r <- build_surface_grid(
    n_x = length(x), window_size = window_size,
    window_increment = window_increment, lag_max = lag_max,
    lag_increment = 1L, lagged = TRUE
  )$n_r
  list(
    observed = obs, null_mean = null_mean, null_sd = null_sd,
    es = es, p = p, n_windows = n_r
  )
}

#' @noRd
.mv_granger_cell <- function(x, y, window_size, window_increment, ar_order, y_surr) {
  surr_res <- wgranger_surrogate(
    x = x, y = y,
    y_surrogates = y_surr,
    window_size = window_size,
    window_increment = window_increment,
    ar_order = ar_order
  )
  # x -> y direction (primary)
  obs_xy <- surr_res$observed_f_xy
  surr_xy <- surr_res$surrogate_f_xy
  null_mean_xy <- mean(surr_xy, na.rm = TRUE)
  null_sd_xy <- stats::sd(surr_xy, na.rm = TRUE)
  es_xy <- if (!is.na(null_sd_xy) && null_sd_xy > 0) {
    (obs_xy - null_mean_xy) / null_sd_xy
  } else {
    NA_real_
  }
  p_xy <- surr_res$p_value_xy

  # y -> x direction
  obs_yx <- surr_res$observed_f_yx
  surr_yx <- surr_res$surrogate_f_yx
  null_mean_yx <- mean(surr_yx, na.rm = TRUE)
  null_sd_yx <- stats::sd(surr_yx, na.rm = TRUE)
  es_yx <- if (!is.na(null_sd_yx) && null_sd_yx > 0) {
    (obs_yx - null_mean_yx) / null_sd_yx
  } else {
    NA_real_
  }
  p_yx <- surr_res$p_value_yx

  n_r <- build_surface_grid(
    n_x = length(x), window_size = window_size,
    window_increment = window_increment, lagged = FALSE
  )$n_r

  list(
    observed = obs_xy, null_mean = null_mean_xy, null_sd = null_sd_xy,
    es = es_xy, p = p_xy,
    observed_yx = obs_yx, null_mean_yx = null_mean_yx, null_sd_yx = null_sd_yx,
    es_yx = es_yx, p_yx = p_yx,
    n_windows = n_r
  )
}


# S3 methods for bsync_multiverse ----------------------------------------------

#' Print method for bsync_multiverse objects
#'
#' @param x A `bsync_multiverse` object.
#' @param ... Additional arguments (not used).
#' @return Returns `x` invisibly.
#' @export
print.bsync_multiverse <- function(x, ...) {
  s <- x$settings
  rb <- x$robustness

  cli::cli_h1("Synchrony Multiverse Analysis ({s$estimator})")
  cli::cli_dl(c(
    "Cells evaluated" = "{rb$n_cells}",
    "Surrogates per cell" = "{s$n_surrogates}",
    "Cells significant (p < .05)" = "{rb$n_significant} ({round(rb$pct_significant * 100, 1)}%)",
    "Median ES" = "{round(rb$median_es, 3)} [IQR: {round(rb$iqr_es, 3)}]",
    "Sign-consistent (sig. cells)" = "{round(rb$sign_consistent * 100, 1)}%"
  ))

  invisible(x)
}

#' Summary method for bsync_multiverse objects
#'
#' @param object A `bsync_multiverse` object.
#' @param ... Additional arguments (not used).
#' @return Returns `object` invisibly.
#' @export
summary.bsync_multiverse <- function(object, ...) {
  print(object)

  cli::cli_h2("Specification Grid")
  cli::cli_text(
    "{nrow(object$grid)} total cells (including {sum(is.na(object$grid$es))} skipped/NA)"
  )

  es_range <- range(object$grid$es, na.rm = TRUE)
  cli::cli_dl(c(
    "ES range"     = "[{round(es_range[1], 3)}, {round(es_range[2], 3)}]",
    "Sample rate"  = "{object$settings$sample_rate} Hz"
  ))

  invisible(object)
}

# Shared windowed-surface infrastructure (M5) ---------------------------------
#
# Centralizes three things that were previously copy-pasted across wcc/wdtw/
# wgranger and their surrogate wrappers:
#
#   1. validate_series()      — x/y/time type + length checks
#   2. validate_window_params() — integerish + positivity checks
#   3. build_surface_grid()   — sole source of the n_r math, w_max = window_size-1,
#                               lag sequence, and short-series abort (Invariant 4)
#   4. new_bsync_surface()    — superclass constructor for all three result objects
#
# The surrogate engine (run_surrogate_engine()) lives in surrogate_engine.R.


# Validation helpers ----------------------------------------------------------

#' @noRd
validate_series <- function(x, y, time = NULL) {
  if (!is.numeric(x)) cli::cli_abort("{.arg x} must be a numeric vector.")
  if (!is.numeric(y)) cli::cli_abort("{.arg y} must be a numeric vector.")
  if (length(x) != length(y)) {
    cli::cli_abort("{.arg x} and {.arg y} must be the same length.")
  }
  if (!is.null(time)) {
    if (!is.numeric(time)) {
      cli::cli_abort("{.arg time} must be a numeric vector.")
    }
    if (length(time) != length(x)) {
      cli::cli_abort("{.arg time} must be the same length as {.arg x}.")
    }
  }
  invisible(NULL)
}

#' @noRd
validate_window_params <- function(window_size, window_increment = 1L,
                                   lag_max = NULL, lag_increment = NULL,
                                   ar_order = NULL, na.rm = NULL) {
  if (!rlang::is_integerish(window_size, n = 1) || window_size <= 0) {
    cli::cli_abort("{.arg window_size} must be a single positive integer.")
  }
  if (!rlang::is_integerish(window_increment, n = 1) || window_increment <= 0) {
    cli::cli_abort("{.arg window_increment} must be a single positive integer.")
  }
  if (!is.null(lag_max)) {
    if (!rlang::is_integerish(lag_max, n = 1) || lag_max <= 0) {
      cli::cli_abort("{.arg lag_max} must be a single positive integer.")
    }
  }
  if (!is.null(lag_increment)) {
    if (!rlang::is_integerish(lag_increment, n = 1) || lag_increment <= 0) {
      cli::cli_abort("{.arg lag_increment} must be a single positive integer.")
    }
  }
  if (!is.null(ar_order)) {
    if (!rlang::is_integerish(ar_order, n = 1) || ar_order <= 0) {
      cli::cli_abort("{.arg ar_order} must be a single positive integer.")
    }
  }
  if (!is.null(na.rm)) {
    if (!rlang::is_logical(na.rm, n = 1)) {
      cli::cli_abort("{.arg na.rm} must be a single logical value.")
    }
  }
  invisible(NULL)
}


# Grid builder ----------------------------------------------------------------

#' Build a windowed-surface (i, tau) grid
#'
#' Sole source of the `n_r` math, the `w_max = window_size - 1` boundary
#' (Invariant 4), the lag sequence, and the short-series abort for all three
#' estimators and their surrogate wrappers. Pass `lagged = FALSE` for Granger
#' (no `tau` dimension).
#'
#' @param n_x Length of the input series.
#' @param window_size Positive integer; number of samples per window.
#' @param window_increment Positive integer; step between window starts.
#' @param lag_max Positive integer; maximum lag (ignored when `lagged = FALSE`).
#' @param lag_increment Positive integer; step between lags (ignored when
#'   `lagged = FALSE`).
#' @param lagged Logical; `TRUE` for WCC/WDTW (has a tau dimension),
#'   `FALSE` for Granger.
#' @return A list with:
#'   - `i_vals`    integer vector of window-start positions (1-based)
#'   - `tau_vals`  integer vector of lags (same length as `i_vals`); `NULL`
#'                 when `lagged = FALSE`
#'   - `w_max`     integer; `window_size - 1` (the C++ boundary value)
#'   - `n_r`       integer; number of distinct window positions
#'   - `lags`      integer vector of unique lags; `NULL` when `lagged = FALSE`
#' @noRd
build_surface_grid <- function(n_x, window_size, window_increment = 1L,
                               lag_max = NULL, lag_increment = 1L,
                               lagged = TRUE) {
  w_max <- window_size - 1L
  w_inc <- window_increment

  if (lagged) {
    tau_max <- lag_max
    tau_inc <- lag_increment
    n_r <- floor((n_x - w_max - 2L * tau_max) / w_inc)

    if (n_r < 1L) {
      cli::cli_abort(c(
        "Series is too short for the requested {.arg window_size} and {.arg lag_max}.",
        "i" = "Need at least {w_max + 1L + 2L * tau_max + 1L} samples; got {n_x}."
      ))
    }

    lags <- seq.int(-tau_max, tau_max, by = tau_inc)
    grid <- base::expand.grid(row = seq_len(n_r), col = seq_along(lags))
    i_vals <- 1L + tau_max + (grid$row - 1L) * w_inc
    tau_vals <- lags[grid$col]

    list(
      i_vals    = i_vals,
      tau_vals  = tau_vals,
      w_max     = w_max,
      n_r       = n_r,
      lags      = lags
    )
  } else {
    n_r <- floor((n_x - w_max) / w_inc)

    if (n_r < 1L) {
      cli::cli_abort(c(
        "Series is too short for the requested {.arg window_size}.",
        "i" = "Need at least {w_max + 1L + 1L} samples; got {n_x}."
      ))
    }

    i_vals <- 1L + (seq_len(n_r) - 1L) * w_inc

    list(
      i_vals   = i_vals,
      tau_vals = NULL,
      w_max    = w_max,
      n_r      = n_r,
      lags     = NULL
    )
  }
}


# Superclass constructor ------------------------------------------------------

#' Construct a bsync_surface result object
#'
#' Adds `"bsync_surface"` to the class vector so shared methods (tidy, glance,
#' as_tibble) can dispatch on it. The estimator-specific leaf class is prepended
#' by the caller.
#'
#' @param x A list with at minimum `results_df`, `settings`, and `aggregate`.
#' @param leaf_class Character; the specific class name, e.g. `"wcc_res"`.
#' @return `x` with class `c(leaf_class, "bsync_surface", "list")`.
#' @noRd
new_bsync_surface <- function(x, leaf_class) {
  stopifnot(is.list(x))
  structure(x, class = c(leaf_class, "bsync_surface", "list"))
}

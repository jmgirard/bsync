# Main Functions ----------------------------------------------------------

#' Windowed Phase Synchrony
#'
#' Conduct a windowed phase-synchrony analysis. Instantaneous phase is
#' extracted once per series with the analytic signal
#' (\code{gsignal::hilbert()}); each window-by-lag cell then summarizes the
#' phase differences into a phase-locking value (PLV) and a mean relative
#' phase.
#'
#' @details
#' The PLV is the modulus of the mean unit phasor of the phase differences —
#' the phase-locking value of Lachaux, Rodriguez, Martinerie, and Varela
#' (1999). Lachaux et al. average across trials at a fixed latency; `wphase()`
#' is a single-trial continuous-data adaptation that averages across the time
#' samples within each sliding window (at each lag offset). PLV is 1 when the
#' phase difference is constant across the window and near 0 when it is
#' uniformly scattered. `rel_phase` is the circular mean of the phase
#' differences (radians); positive values mean `x`'s phase is ahead of `y`'s.
#'
#' **Narrowband caveat (stated loudly, per the source method):** instantaneous
#' phase is well-defined for narrowband signals. Lachaux et al. band-pass
#' filter before phase extraction; `wphase()` leaves filtering to you (e.g.
#' [smooth_signal()] or an external band-pass) rather than silently changing
#' the basis. Beware that this failure mode is **anti-conservative**:
#' broadband or low-frequency-dominated input *inflates* PLV toward 1 (two
#' independent AR(0.98) series measure a mean PLV near 0.6 at
#' `window_size = 96`, against roughly 0.12 for white noise), because slowly
#' drifting phases stay aligned across a short window. High PLV on
#' unfiltered signals is not evidence of synchrony; band-pass first and use
#' [wphase_surrogate()] for significance.
#'
#' **NA policy:** the analytic signal is FFT-based, so a single `NA` corrupts
#' every phase estimate; unlike [wcc()], there is no coherent per-window
#' `na.rm`. `wphase()` therefore aborts on `NA`-containing input — impute with
#' [impute_ts_gaps()] or trim with [trim_edges()] (`pad_na = FALSE`) first.
#'
#' @param x A numeric vector containing a time series (same length as `y`).
#' @param y A numeric vector containing a time series (same length as `x`).
#' @param time An optional numeric vector of timestamps for the data, same
#'   length as `x` and `y`; if provided, window indices in the results are
#'   mapped to these timestamps (and it must not contain missing values).
#'   Default is `NULL`.
#' @param window_size A positive integer indicating the size of each window in
#'   samples.
#' @param lag_max A positive integer indicating the maximum lag (in samples)
#'   to try between the `x` and `y` windows.
#' @param window_increment A positive integer indicating the number of samples
#'   between successive windows. (default = `1`)
#' @param lag_increment A positive integer indicating the number of samples
#'   between successive lags. (default = `1`)
#' @param statistic A character string naming how to aggregate the PLV
#'   surface into a single number. `"mean_plv"` (default) takes the mean PLV
#'   over **all** windows and lags. It summarizes the whole surface and is the
#'   value `wphase()` reported before this argument existed. `"peak"` takes
#'   the largest PLV across lags **within each window**, then averages those
#'   per-window values. This is the rMEA *best-lag* convention (Boker et al.,
#'   2002) that [wcc()] uses for its own `"peak"` statistic, applied to PLV.
#'   It suits a dyad whose lead–lag shifts across windows, where the mean
#'   over all lags dilutes the lag at which phases lock. Pass the same value
#'   to [wphase_surrogate()] so that the null distribution matches.
#' @return A list object of class `"wphase_res"` (a `bsync_surface`) with
#'   `results_df` (columns `i`, `tau`, `plv`, `rel_phase`), `aggregate`
#'   (named `mean_plv` or `peak`, after `statistic`), and `settings`.
#' @seealso [wphase_surrogate()] for the matched surrogate test, [pick_optima()]
#'   for per-window optima, and `vignette("wphase-workflow")` for a complete
#'   workflow with band-pass filtering.
#' @examples
#' # Windowed phase synchrony on the bundled simulated dyad
#' wphase_res <- wphase(
#'   x = sim_dyad$x_A,
#'   y = sim_dyad$x_B,
#'   window_size = 96,
#'   lag_max = 10
#' )
#' wphase_res
#'
#' # Average the best-lag PLV of each window instead
#' wphase(sim_dyad$x_A, sim_dyad$x_B, window_size = 96, lag_max = 10, statistic = "peak")
#' @md
#' @export
wphase <- function(
  x,
  y,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1,
  statistic = c("mean_plv", "peak")
) {
  statistic <- match.arg(statistic)
  validate_series(x, y, time)
  validate_window_params(
    window_size, window_increment,
    lag_max = lag_max, lag_increment = lag_increment
  )
  if (anyNA(x) || anyNA(y) || (!is.null(time) && anyNA(time))) {
    cli::cli_abort(c(
      "{.arg x} and {.arg y} must not contain missing values.",
      "i" = "The analytic signal (FFT-based Hilbert transform) is corrupted by any NA.",
      "i" = "Impute with {.fn impute_ts_gaps} or trim with {.fn trim_edges} (pad_na = FALSE) first."
    ))
  }

  x <- as.double(x)
  y <- as.double(y)

  settings <- list(
    window_size = window_size,
    window_increment = window_increment,
    lag_max = lag_max,
    lag_increment = lag_increment,
    statistic = statistic,
    has_time = !is.null(time)
  )

  surface <- create_wphase_df(
    x = x,
    y = y,
    time = time,
    settings = settings
  )
  results_df <- surface$results_df

  agg_val <- wphase_aggregate(
    plv = results_df$plv,
    window_id = surface$window_id,
    statistic = statistic
  )

  out <- list(
    results_df = results_df,
    aggregate  = stats::setNames(agg_val, statistic),
    settings   = settings
  )

  new_bsync_surface(out, "wphase_res")
}

# S3 Methods --------------------------------------------------------------

#' Print method for wphase_res objects
#'
#' @param x An object of class "wphase_res".
#' @param ... Additional arguments (not used).
#' @return Returns `x` invisibly.
#' @export
print.wphase_res <- function(x, ...) {
  s <- x$settings
  n_windows <- length(unique(x$results_df$i))
  n_lags <- length(unique(x$results_df$tau))

  cli::cli_h1("Windowed Phase Synchrony Analysis")

  cli::cli_dl(stats::setNames(
    c(
      "{n_windows}",
      "{n_lags}",
      "{s$window_size}",
      "{s$lag_max}",
      "{round(x$aggregate[[1]], 4)}"
    ),
    c(
      "Total Windows",
      "Total Lags Tested",
      "Window Size",
      "Max Lag",
      wphase_agg_label(s$statistic)
    )
  ))

  invisible(x)
}

#' Summary method for wphase_res objects
#'
#' @param object An object of class "wphase_res".
#' @param ... Additional arguments (not used).
#' @return Returns `object` invisibly.
#' @export
summary.wphase_res <- function(object, ...) {
  print(object)

  cli::cli_h2("Phase-Locking Value Distribution")
  plv_vals <- object$results_df$plv
  q_vals <- stats::quantile(
    plv_vals,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    na.rm = TRUE
  )

  print(round(q_vals, 4))

  n_na <- sum(is.na(plv_vals))
  if (n_na > 0) {
    cli::cli_alert_warning("{n_na} missing value{?s} (NA) detected.")
  }

  invisible(object)
}

# Internal Helpers --------------------------------------------------------

#' @noRd
extract_phase <- function(v) {
  Arg(gsignal::hilbert(v))
}

#' @noRd
create_wphase_df <- function(x, y, time = NULL, settings) {
  grid <- build_surface_grid(
    n_x = length(x),
    window_size = settings$window_size,
    window_increment = settings$window_increment,
    lag_max = settings$lag_max,
    lag_increment = settings$lag_increment,
    lagged = TRUE
  )

  core <- calc_wphase_cpp(
    phi_x    = extract_phase(x),
    phi_y    = extract_phase(y),
    i_vals   = grid$i_vals,
    tau_vals = grid$tau_vals,
    w_max    = grid$w_max
  )

  results_df <- data.frame(
    i = grid$i_vals,
    tau = grid$tau_vals,
    plv = core$plv,
    rel_phase = core$rel_phase
  )

  if (!is.null(time)) {
    results_df$i <- time[results_df$i]
  }

  # window_id holds the grid positions from before the time mapping, so that
  # tied timestamps cannot merge two windows (wphase_surrogate() groups by the
  # same grid positions).
  list(results_df = results_df, window_id = grid$i_vals)
}

# Console label of the aggregate, shared by the wphase_res and wphase_surr
# print methods.
#' @noRd
wphase_agg_label <- function(statistic) {
  if (identical(statistic, "peak")) "Mean Peak PLV" else "Mean PLV"
}

# The one aggregate for the observed surface and every surrogate (Invariant
# 2). "mean_plv": mean PLV over all windows and lags. "peak": per-window
# maximum PLV across lags, averaged over windows -- the wcc "peak" definition
# (wcc_aggregate()) applied to PLV. window_id: grid window positions.
#' @noRd
wphase_aggregate <- function(plv, window_id, statistic = "mean_plv") {
  if (statistic == "mean_plv") {
    base::mean(plv, na.rm = TRUE)
  } else {
    peaks <- tapply(plv, window_id, function(v) {
      if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    })
    base::mean(peaks, na.rm = TRUE)
  }
}

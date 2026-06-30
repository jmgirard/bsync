# Auto-tune WCC parameters (M6) -----------------------------------------------
#
# autotune_wcc() is a thin wrapper over synchrony_multiverse() that applies a
# gated stability-penalized selection rule across a dyad_list:
#
#   GATE: specification must be significant in >= sig_pct of dyads
#   SCORE: median ES across dyads - iqr_penalty * IQR(ES across dyads)
#
# select_specification() is the internal helper that runs on a list of
# bsync_multiverse objects (one per dyad) and returns the winning row index.
#
# Invariant 6: stochastics respect set.seed() / future.seed; no internal
# reseeding. The per-dyad synchrony_multiverse() calls use future.apply
# internally, which honors future.seed.


#' Auto-Tune WCC Parameters for a Multi-Dyad Dataset
#'
#' Selects Windowed Cross-Correlation hyperparameters that are both detectable
#' (significant vs. the null) and stable (consistent) across a collection of
#' dyads. Internally calls [synchrony_multiverse()] on each dyad and applies
#' a gated stability-penalized selection rule via [select_specification()].
#'
#' @details
#' **Why cross-dyad stability?** A parameter set that maximizes raw synchrony
#' for one dyad may simply match that dyad's autocorrelation structure. The
#' matched-null surrogate controls for autocorrelation within a dyad (Invariant
#' 2), but the *best* parameters should also replicate across dyads with
#' structurally different signals -- hence the multi-dyad stability criterion.
#'
#' **Selection rule.** Cells pass a detectability gate (significant in at least
#' `sig_pct` of dyads). Among passing cells, the score is
#' `median(ES) - iqr_penalty * IQR(ES)` across dyads, penalizing spread.
#' If no cell passes the gate, a warning is issued and the highest-median-ES
#' cell is returned (soft fallback).
#'
#' **Dyad sampling.** If `length(dyad_list) > n_tune_dyads`, a random sample
#' of `n_tune_dyads` dyads is used for speed; call `set.seed()` beforehand for
#' reproducibility.
#'
#' @param dyad_list A list of data frames or lists. Each element represents one
#'   dyad and must have at least two numeric columns (or two named list elements
#'   `x` and `y`) containing the two time series.
#' @param sample_rate Single positive number; sampling rate in Hz, used to
#'   convert `window_sec` and `lag_sec` to samples.
#' @param window_sec Numeric vector; window size(s) in seconds to sweep.
#'   Use [suggest_wcc_params()] on a representative dyad to find a principled
#'   starting range.
#' @param lag_sec Numeric vector; max lag(s) in seconds. Default `NULL` uses
#'   `window_sec / 2` per cell (the SUSY reliability ceiling).
#' @param increment_pct Numeric; window increment as a fraction of window size
#'   (e.g., `0.1` = 10\% step). Default is `0.1`.
#' @param statistic Character; WCC aggregate statistic. Default `"mean_abs_z"`.
#' @param surrogate_method Character; surrogate generator: `"phase"` (default)
#'   or `"circular"`.
#' @param n_surrogates Single positive integer; surrogates per cell per dyad.
#'   Default `100`. Increase to >= 1000 for reporting.
#' @param n_tune_dyads Maximum number of dyads to use. If
#'   `length(dyad_list) > n_tune_dyads`, a random sample is taken. Default
#'   `30`.
#' @param sig_pct Detectability gate: minimum proportion of dyads in which a
#'   cell must be significant (p < .05). Default `0.5`.
#' @param iqr_penalty Penalty weight on cross-dyad IQR of ES. Score =
#'   `median(ES) - iqr_penalty * IQR(ES)`. Default `0.5`.
#' @return An object of class `bsync_autotune` (a named list with a tidy
#'   [print()][print.bsync_autotune()] method) containing:
#'   \describe{
#'     \item{`window_size`}{Selected window size in samples.}
#'     \item{`lag_max`}{Selected max lag in samples.}
#'     \item{`window_increment`}{Selected window increment in samples.}
#'     \item{`lag_increment`}{`1L` (standard lag increment).}
#'     \item{`window_sec`}{Selected window size in seconds.}
#'     \item{`lag_sec`}{Selected max lag in seconds.}
#'     \item{`sig_rate`}{Proportion of dyads where selected cell was significant.}
#'     \item{`median_es`}{Median ES across dyads for the selected cell.}
#'     \item{`iqr_es`}{IQR of ES across dyads for the selected cell.}
#'     \item{`score`}{Selection score for the chosen cell.}
#'     \item{`n_dyads`}{Number of dyads used for tuning.}
#'     \item{`n_cells_gated`}{Number of cells that passed the detectability gate.}
#'     \item{`dyad_multiverses`}{List of `bsync_multiverse` objects, one per dyad.}
#'   }
#' @seealso [synchrony_multiverse()], [suggest_wcc_params()],
#'   [select_specification()]
#' @examples
#' \donttest{
#' # Tune across a small multi-dyad list (here three copies of one dyad).
#' # Small surrogate count for a fast example; use >= 1000 for reporting.
#' dyads <- replicate(
#'   3,
#'   list(x = sim_dyad$x_A, y = sim_dyad$x_B),
#'   simplify = FALSE
#' )
#' tuned <- autotune_wcc(
#'   dyad_list = dyads,
#'   sample_rate = 80,
#'   window_sec = c(1, 2, 4),
#'   lag_sec = 1,
#'   n_surrogates = 30
#' )
#' tuned
#' }
#' @export
autotune_wcc <- function(
  dyad_list,
  sample_rate,
  window_sec,
  lag_sec = NULL,
  increment_pct = 0.1,
  statistic = "mean_abs_z",
  surrogate_method = "phase",
  n_surrogates = 100L,
  n_tune_dyads = 30L,
  sig_pct = 0.5,
  iqr_penalty = 0.5
) {
  if (!is.list(dyad_list) || length(dyad_list) < 1) {
    cli::cli_abort("{.arg dyad_list} must be a non-empty list.")
  }
  if (!is.numeric(sample_rate) || length(sample_rate) != 1 || sample_rate <= 0) {
    cli::cli_abort("{.arg sample_rate} must be a single positive number.")
  }
  if (!is.numeric(window_sec) || length(window_sec) < 1 || any(window_sec <= 0)) {
    cli::cli_abort("{.arg window_sec} must be a positive numeric vector.")
  }
  if (!rlang::is_integerish(n_tune_dyads, n = 1) || n_tune_dyads < 1) {
    cli::cli_abort("{.arg n_tune_dyads} must be a single positive integer.")
  }
  if (!is.numeric(sig_pct) || length(sig_pct) != 1 ||
    sig_pct < 0 || sig_pct > 1) {
    cli::cli_abort("{.arg sig_pct} must be a single number in [0, 1].")
  }
  if (!is.numeric(iqr_penalty) || length(iqr_penalty) != 1 || iqr_penalty < 0) {
    cli::cli_abort("{.arg iqr_penalty} must be a single non-negative number.")
  }

  # Default lag_sec: SUSY ceiling (window / 2)
  lag_sec_use <- lag_sec %||% (window_sec / 2)

  # Sample dyads
  n_total <- length(dyad_list)
  n_use <- min(n_tune_dyads, n_total)
  if (n_use < n_total) {
    tune_idx <- sample.int(n_total, n_use)
    cli::cli_inform("Sampling {n_use} of {n_total} dyads for tuning.")
  } else {
    tune_idx <- seq_len(n_total)
  }
  tune_list <- dyad_list[tune_idx]

  cli::cli_inform("Running synchrony_multiverse() on {n_use} dyad(s) \\
    ({length(window_sec)} window x {length(lag_sec_use)} lag cells each)...")

  # Run multiverse on each dyad
  mv_list <- lapply(tune_list, function(dyad) {
    xy <- .extract_xy(dyad)
    synchrony_multiverse(
      x = xy$x,
      y = xy$y,
      estimator = "wcc",
      sample_rate = sample_rate,
      window_sec = window_sec,
      lag_sec = lag_sec_use,
      increment_pct = increment_pct,
      statistic = statistic,
      surrogate_method = surrogate_method,
      n_surrogates = n_surrogates
    )
  })

  # Apply selection rule
  sel <- select_specification(
    mv_list,
    sig_pct = sig_pct, iqr_penalty = iqr_penalty
  )
  best <- sel$best_row

  # Assemble result
  result <- list(
    window_size      = best$window_size,
    lag_max          = best$lag_max,
    window_increment = best$window_increment,
    lag_increment    = 1L,
    window_sec       = best$window_sec,
    lag_sec          = best$lag_sec,
    sig_rate         = sel$sig_rate,
    median_es        = sel$median_es,
    iqr_es           = sel$iqr_es,
    score            = sel$score,
    n_dyads          = n_use,
    n_cells_gated    = sel$n_gated,
    dyad_multiverses = mv_list
  )
  class(result) <- "bsync_autotune"

  result
}


#' Print method for bsync_autotune objects
#'
#' @param x An object of class `bsync_autotune` (from [autotune_wcc()]).
#' @param ... Additional arguments (not used).
#' @return Returns `x` invisibly.
#' @export
print.bsync_autotune <- function(x, ...) {
  cli::cli_h2("Auto-Tune Result")
  cli::cli_dl(c(
    "Window size" = "{x$window_size} samples ({round(x$window_sec, 2)} s)",
    "Max lag"     = "{x$lag_max} samples ({round(x$lag_sec, 2)} s)",
    "Increment"   = "{x$window_increment} samples",
    "Sig. rate"   = "{round(x$sig_rate * 100, 1)}% of dyads",
    "Median ES"   = "{round(x$median_es, 3)} (IQR = {round(x$iqr_es, 3)})"
  ))
  cli::cli_alert_info(
    "Tuned over {x$n_dyads} dyad{?s}; {x$n_cells_gated} cell{?s} passed the \\
     detectability gate. Per-dyad multiverses in {.field $dyad_multiverses}."
  )
  invisible(x)
}


#' Select the Best Specification from a Multi-Dyad Multiverse
#'
#' Applies the gated stability-penalized selection rule to a list of
#' `bsync_multiverse` objects (one per dyad) and returns the winning cell.
#'
#' @param mv_list A list of `bsync_multiverse` objects, all run with the same
#'   parameter grid (i.e., all produced by [synchrony_multiverse()] with
#'   identical `window_sec`, `lag_sec`, `increment_pct`, `statistic`, and
#'   `surrogate_method` arguments).
#' @param sig_pct Minimum proportion of dyads in which a cell must be
#'   significant (p < .05) to pass the detectability gate. Default `0.5`.
#' @param iqr_penalty Penalty weight on cross-dyad IQR of ES in the score
#'   `median(ES) - iqr_penalty * IQR(ES)`. Default `0.5`.
#' @return A list with `best_row` (one-row tibble from the grid), `sig_rate`,
#'   `median_es`, `iqr_es`, `score`, and `n_gated` for the selected cell.
#' @seealso [autotune_wcc()], [synchrony_multiverse()]
#' @examples
#' \donttest{
#' # Build one multiverse per dyad, then pick the most robust specification.
#' # Small surrogate count for a fast example; use >= 1000 for reporting.
#' mv_list <- lapply(seq_len(3), function(i) {
#'   synchrony_multiverse(
#'     x = sim_dyad$x_A,
#'     y = sim_dyad$x_B,
#'     estimator = "wcc",
#'     sample_rate = 80,
#'     window_sec = c(1, 2, 4),
#'     lag_sec = 1,
#'     n_surrogates = 30
#'   )
#' })
#' select_specification(mv_list)
#' }
#' @export
select_specification <- function(mv_list, sig_pct = 0.5, iqr_penalty = 0.5) {
  if (!is.list(mv_list) || length(mv_list) < 1 ||
    !all(vapply(mv_list, inherits, logical(1), "bsync_multiverse"))) {
    cli::cli_abort(
      "{.arg mv_list} must be a non-empty list of {.cls bsync_multiverse} objects."
    )
  }

  grids <- lapply(mv_list, function(mv) mv$grid)
  n_dyads <- length(grids)
  n_cells <- nrow(grids[[1]])

  sig_rate <- numeric(n_cells)
  median_es <- numeric(n_cells)
  iqr_es <- numeric(n_cells)

  for (j in seq_len(n_cells)) {
    es_j <- vapply(grids, function(g) g$es[j], numeric(1))
    p_j <- vapply(grids, function(g) g$p[j], numeric(1))
    ok <- !is.na(es_j) & !is.na(p_j)
    if (!any(ok)) {
      sig_rate[j] <- NA_real_
      median_es[j] <- NA_real_
      iqr_es[j] <- NA_real_
    } else {
      sig_rate[j] <- mean(p_j[ok] < 0.05)
      median_es[j] <- stats::median(es_j[ok])
      iqr_es[j] <- stats::IQR(es_j[ok])
    }
  }

  gated <- !is.na(sig_rate) & sig_rate >= sig_pct

  if (!any(gated)) {
    cli::cli_warn(
      "No specification passed the detectability gate ({.val {sig_pct}} significance rate). \\
      Falling back to highest median ES."
    )
    gated <- !is.na(median_es)
  }

  score <- median_es - iqr_penalty * iqr_es
  score[!gated] <- NA_real_

  best_idx <- which.max(score)

  list(
    best_row  = grids[[1]][best_idx, ],
    sig_rate  = sig_rate[best_idx],
    median_es = median_es[best_idx],
    iqr_es    = iqr_es[best_idx],
    score     = score[best_idx],
    n_gated   = sum(gated, na.rm = TRUE)
  )
}


# .extract_xy() ----------------------------------------------------------------
# Internal helper: pull x and y from a dyad element (data.frame with 2+ cols
# or list with $x and $y).

.extract_xy <- function(dyad) {
  if (is.data.frame(dyad)) {
    if (ncol(dyad) < 2) {
      cli::cli_abort("Each dyad data frame must have at least two columns.")
    }
    list(x = dyad[[1]], y = dyad[[2]])
  } else if (is.list(dyad)) {
    if (!is.null(dyad$x) && !is.null(dyad$y)) {
      list(x = dyad$x, y = dyad$y)
    } else if (length(dyad) >= 2) {
      list(x = dyad[[1]], y = dyad[[2]])
    } else {
      cli::cli_abort("Each dyad must be a data frame or a list with at least two elements.")
    }
  } else {
    cli::cli_abort("Each dyad must be a data frame or a list.")
  }
}

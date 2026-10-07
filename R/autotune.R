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
#' @section Avoiding circular inference:
#' `autotune_wcc()` **selects** the parameters that maximize the stability-
#' penalized effect size. If you then run a confirmatory synchrony test with
#' those parameters **on the same dyads**, the reported effect size and p-value
#' are inflated -- the parameters were chosen to maximize exactly that quantity
#' (a form of double-dipping / circular analysis). Treat the tuned parameters as
#' a *hypothesis*, not a validated result. To draw defensible conclusions,
#' choose one of:
#' \itemize{
#'   \item **Split the sample.** Tune on a pilot/exploratory subset of dyads and
#'     confirm on a disjoint held-out subset that was not passed to
#'     `autotune_wcc()`.
#'   \item **Report the multiverse, not the winner.** Instead of a single tuned
#'     specification, report the robustness of synchrony across the whole grid
#'     via [synchrony_multiverse()] (`pct_significant`, `sign_consistent`, the ES
#'     distribution). A result that holds only at the single best-tuned cell is
#'     weak evidence; one that holds across the grid is strong.
#' }
#' The per-dyad multiverses are returned in `$dyad_multiverses` precisely so you
#' can inspect that full distribution rather than collapsing to the point
#' estimate.
#'
#' @param dyad_list A list of data frames or lists, one element per dyad. If
#'   an element's names contain \code{x} and \code{y} exactly once each, those
#'   two elements are read by name, in any position. If its names contain
#'   neither \code{x} nor \code{y}, or it has no names, its first two elements
#'   are read as \code{x} and \code{y}. Any other use of the names \code{x}
#'   and \code{y} is an error. Names are never matched partially. The names
#'   and shape of every dyad are checked before any dyads are sampled.
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
#' @param surrogate_method Character; surrogate generator: \code{"phase"}
#'   (default), \code{"circular"}, or \code{"iaaft"}. \code{"iaaft"} needs
#'   every \code{y} without missing values. See
#'   \code{\link{synchrony_multiverse}}.
#' @param n_surrogates Single positive integer; surrogates per cell per dyad.
#'   Default `100`. Increase to >= 1000 for reporting. Below 19 the call
#'   warns once, because no cell can then reach p <= .05, so the
#'   detectability gate passes only if `sig_pct = 0`.
#' @param n_tune_dyads Maximum number of dyads to use. If
#'   `length(dyad_list) > n_tune_dyads`, a random sample is taken. Default
#'   `30`.
#' @param sig_pct Detectability gate: minimum proportion of dyads in which a
#'   cell must be significant (p <= .05). Default `0.5`.
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
#'   [select_specification()]; \code{\link{generate_surrogate_iaaft}} for
#'   the \code{"iaaft"} method
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
  if (is.data.frame(dyad_list)) {
    cli::cli_abort(
      "{.arg dyad_list} must be a list of dyads, not a single data frame. \\
      Wrap one data frame per dyad in a list."
    )
  }
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
  if (!rlang::is_integerish(n_surrogates, n = 1) || is.na(n_surrogates) ||
    n_surrogates < 1) {
    cli::cli_abort("{.arg n_surrogates} must be a single positive integer.")
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

  # Read every dyad before sampling, so a malformed dyad aborts the call even
  # when the sample would leave it out, and the abort names its true index.
  xy_list <- lapply(
    seq_along(dyad_list),
    function(i) .extract_xy(dyad_list[[i]], i)
  )
  names(xy_list) <- names(dyad_list)

  # Sample dyads
  n_total <- length(dyad_list)
  n_use <- min(n_tune_dyads, n_total)
  if (n_use < n_total) {
    tune_idx <- sample.int(n_total, n_use)
    cli::cli_inform("Sampling {n_use} of {n_total} dyads for tuning.")
  } else {
    tune_idx <- seq_len(n_total)
  }
  tune_list <- xy_list[tune_idx]

  # IAAFT cannot take NA, so name every sampled dyad that has one up front.
  if ("iaaft" %in% surrogate_method) {
    has_na <- vapply(tune_list, function(xy) anyNA(xy$y), logical(1))
    if (any(has_na)) {
      bad <- tune_idx[has_na]
      n_bad <- length(bad)
      cli::cli_abort(c(
        "{.code surrogate_method = \"iaaft\"} needs every {.arg y} without \\
        missing values.",
        "x" = "{cli::qty(n_bad)}Dyad{?s} {bad} of {.arg dyad_list} \\
        {cli::qty(n_bad)}{?has/have} missing values in {.arg y}.",
        "i" = "Fill gaps first, for example with {.fn impute_ts_gaps}."
      ))
    }
  }

  cli::cli_inform("Running synchrony_multiverse() on {n_use} dyad(s) \\
    ({length(window_sec)} window x {length(lag_sec_use)} lag cells each)...")

  # Give the few-surrogates warning once here, and muffle the copy that
  # each per-dyad synchrony_multiverse() call would give.
  warn_few_surrogates(n_surrogates)

  # Run multiverse on each dyad. IAAFT non-convergence is counted per dyad
  # and reported once after the loop.
  n_unconverged_dyads <- 0L
  mv_list <- lapply(tune_list, function(xy) {
    unconverged <- FALSE
    res <- withCallingHandlers(
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
      ),
      bsync_few_surrogates = function(w) invokeRestart("muffleWarning"),
      bsync_iaaft_unconverged = function(w) {
        unconverged <<- TRUE
        invokeRestart("muffleWarning")
      }
    )
    if (unconverged) n_unconverged_dyads <<- n_unconverged_dyads + 1L
    res
  })
  if (n_unconverged_dyads > 0) {
    cli::cli_warn(
      c(
        "Some IAAFT surrogates did not converge in {n_unconverged_dyads} of \\
        {n_use} dyad{?s}.",
        "i" = "They still match the spectrum of {.arg y} exactly. Their \\
        values are a little less close to the values of {.arg y}."
      ),
      class = "bsync_iaaft_unconverged"
    )
  }

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
  cli::cli_alert_warning(
    "These parameters were {.emph selected} to maximize effect size. Re-testing \\
     synchrony with them on the same dyads is circular; confirm on held-out \\
     dyads or report the multiverse robustness. See {.help autotune_wcc}."
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
#'   significant (p <= .05) to pass the detectability gate. Default `0.5`.
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
      sig_rate[j] <- mean(is_significant(p_j[ok]))
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
# Internal helper: pull x and y from one dyad of a dyad_list (a data frame or
# a list). One rule for both shapes, with exact name matching only:
#   - `x` and `y` each named exactly once -> read by name, in any position;
#   - neither `x` nor `y` named (or no names) -> the first two elements;
#   - any other count of exact `x`/`y` names -> abort (no silent role swap).
# `index` is the dyad's position in the user's dyad_list, named in every
# abort so the user can find the bad dyad.

.extract_xy <- function(dyad, index) {
  if (!is.list(dyad)) {
    cli::cli_abort(c(
      "Dyad {index} in {.arg dyad_list} must be a data frame or a list.",
      "i" = "If you passed a single dyad, wrap it in a list of dyads."
    ))
  }
  nm <- names(dyad)
  n_x <- sum(nm %in% "x")
  n_y <- sum(nm %in% "y")
  if (n_x == 1 && n_y == 1) {
    return(list(
      x = dyad[[which(nm %in% "x")]],
      y = dyad[[which(nm %in% "y")]]
    ))
  }
  if (n_x > 0 || n_y > 0) {
    cli::cli_abort(c(
      "Dyad {index} in {.arg dyad_list} names {.val x} {n_x} time{?s} and \\
      {.val y} {n_y} time{?s}.",
      "i" = "Name each exactly once, or use neither name to read the first \\
      two elements by position."
    ))
  }
  if (length(dyad) < 2) {
    cli::cli_abort(
      "Dyad {index} in {.arg dyad_list} must have at least two elements \\
      (columns)."
    )
  }
  list(x = dyad[[1]], y = dyad[[2]])
}

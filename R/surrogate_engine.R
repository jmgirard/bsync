# Shared surrogate engine (M5) -----------------------------------------------
#
# run_surrogate_engine() is the single implementation of the surrogate loop.
# All three surrogate wrappers (wcc_surrogate, wdtw_surrogate, wgranger_surrogate)
# delegate to it.
#
# Design principles (M6 multiverse efficiency, Invariant 7):
#   - Accepts a prebuilt grid (from build_surface_grid()) and a prebuilt
#     y_surrogates matrix so both can be hoisted out of an outer parameter-cell
#     loop in M6 (surrogates depend only on surrogate_method, not on
#     window_size/lag_max/increment).
#   - aggregate-only path: compute_fn(x, y_col, grid) returns a numeric scalar
#     or named numeric — never a full results_df. The heavy per-cell data.frame
#     assembly (create_*_df) is not called on the surrogate path (Invariant 7).
#   - Tail direction and p-value computation are the caller's responsibility
#     (they need the observed stat, which the engine does not hold).

#' Run the shared surrogate significance engine
#'
#' @param x Numeric vector (the reference series); pre-processed (scaled etc.)
#'   as needed by `compute_fn`.
#' @param y_surrogates Numeric matrix; each column is one surrogate for `y`.
#' @param grid Output of `build_surface_grid()`.
#' @param compute_fn `function(x, y_col, grid)` → numeric scalar or named
#'   numeric. Must not construct a full `results_df`; aggregate only.
#' @param fun_value Template for `future_vapply`'s `FUN.VALUE`: use
#'   `numeric(1)` for scalar aggregates (WCC/WDTW), `numeric(2)` for Granger's
#'   two-direction aggregate.
#' @return Numeric vector (fun_value = numeric(1)) or named matrix
#'   (n_stats × n_surrogates, fun_value = numeric(k>1)). P-value computation is
#'   the caller's responsibility.
#' @noRd
run_surrogate_engine <- function(x, y_surrogates, grid, compute_fn,
                                 fun_value = numeric(1)) {
  future.apply::future_vapply(
    seq_len(ncol(y_surrogates)),
    function(idx) compute_fn(x, y_surrogates[, idx], grid),
    FUN.VALUE = fun_value,
    future.seed = TRUE
  )
}

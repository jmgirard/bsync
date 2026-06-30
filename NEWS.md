# bsync 0.0.0.9000

* M1: `na.rm` is now honored in WCC. `calc_wcc_cpp()` gains an `na_rm` parameter;
  `na.rm = FALSE` in `wcc()` / `wcc_surrogate()` now returns `NA` for any window
  containing a missing value instead of silently using pairwise-complete pairs.

* M1: Window-size semantics fixed to exactly `window_size` samples. `w_max` is now
  set to `window_size - 1` at the C++ boundary in all three estimators and surrogate
  wrappers; previously each window spanned `window_size + 1` samples (off-by-one).

* M1: Short-series robustness. All three `create_*_df()` builders and surrogate
  grid builders now call `cli::cli_abort()` when the series is too short for the chosen
  `window_size` / `lag_max`; index sequences use `seq_len()` throughout.

* M1: `evaluate_signal_power()` no longer auto-prints its plot. The `plot` argument
  has been removed; call `plot()` on the returned `"signal_power_res"` object instead,
  eliminating the `Rplots.pdf` side effect in non-interactive contexts.

* M1: Condition style unified in `R/impute.R` and `R/surrogate_generation.R`. Legacy
  `stop()` / `warning()` calls replaced with `cli::cli_abort()` / `cli::cli_warn()`.

* M1: `leadership_asymmetry()` docs now state the centered sliding-window semantics
  explicitly; `min_valid` is validated as a single positive integer.

* M1: `suggest_wcc_params()` now documents the 4-cycles-per-window heuristic
  (`window_size = round(event_duration_sec * 4 * sample_rate)`) in `@details`.

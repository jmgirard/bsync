# bsync 0.0.0.9000

## M3 — CRAN readiness

* **Build artifacts removed from version control.** `src/*.o`, `src/bsync.{so,dll}`,
  `.DS_Store`, and `tests/testthat/Rplots.pdf` were untracked via `git rm --cached`;
  `.gitignore` now carries `*.o`, `*.so`, `*.dll`, `*.dylib`, and `Rplots.pdf` patterns.

* **Documentation complete.** Added `@return` to all 18 exported `print()`, `plot()`,
  and `summary()` methods (`print.wcc_res`, `summary.wcc_res`, `print.wdtw_res`,
  `print/summary.wgranger_res`, `print/summary.wcc_optima`, `print/summary.wdtw_optima`,
  `print.bsync_lai`, `plot.bsync_lai`, `plot.wcc_res`, `plot.wdtw_res`,
  `plot.wgranger_res`, `plot_optima_overlay`, `print.wcc_surr`, `print.wdtw_surr`,
  `print.wgranger_surr`). `R CMD check --as-cran` now reports 0 errors / 0 warnings
  / 0 notes.

* **DESCRIPTION updated.** Description field refreshed to cover all three windowed
  estimators (WCC, WDTW, Granger), surrogate significance testing (circular-shift and
  phase-randomization generators), optima picking, leadership asymmetry index, and the
  full preprocessing pipeline. `Language: en-US` field added.

* **Spell-check clean.** `inst/WORDLIST` created with domain terms, acronyms, and
  proper names; `spelling::spell_check_package()` returns no errors.

* **README corrected.** The surrogate example chunk now uses the correct two-step API:
  `generate_surrogate_circular()` to build the null matrix, then `wcc_surrogate()` with
  `y_surrogates`. `pkgdown::check_pkgdown()` passes.

## M2 — Efficiency

* **WCC core rewritten to an NA-aware prefix-sum algorithm.** `calc_wcc_cpp()` now
  preprocesses six masked prefix-sum arrays per lag in O(n) and evaluates each window
  in O(1), replacing the prior O(w\_max) inner loop. Both `na.rm` modes are preserved
  exactly. Speedup on typical configurations: 5–25× depending on `window_size` and
  `lag_max`; see `bench/RESULTS.md` for measured timings.

* **OpenMP removed.** The prefix-sum speedup makes OpenMP unnecessary for the WCC core.
  `SHLIB_OPENMP_CXXFLAGS` has been stripped from `src/Makevars` and `src/Makevars.win`;
  the package is serial-by-default and fully reproducible on all platforms.

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

# bsync 0.0.0.9000

## M7 — Documentation & messaging overhaul (Phase A)

* **New `vignette("bsync")`** ("Get started") provides a full end-to-end workflow
  walkthrough using `sim_dyad`, a WCC / WDTW / WGC estimator-choice decision table,
  and a reading map into the six deep-dive articles. This is now the recommended
  entry point for new users.

* **Structured pkgdown site**: articles are now grouped into "Get started", "Core
  estimators", and "Going deeper" sections; the function reference index is organized
  into eight thematic groups (Estimators, Surrogate testing, Optima & leadership,
  Parameter guidance, Preprocessing x2, Tidy interface, Data).

* **`choosing-parameters` vignette reworked** — retitled "Choosing Analysis Parameters"
  (was WCC-only); the `synchrony_multiverse()` / `autotune_wcc()` / `select_specification()`
  examples now run on `sim_dyad` so the Simonsohn specification-curve plot and
  `glance()` robustness summary render; WDTW and Granger `estimator=` support shown.

* **Cross-linking**: `wdtw-workflow`, `wgranger-workflow`, and `determine-downsampling`
  vignettes now carry "See also / Next steps" footers consistent with `wcc-workflow`.

* **README refreshed**: Granger causality added to the headline scope sentence;
  parameter-guidance tools (`suggest_wcc_params()`, `synchrony_multiverse()`,
  `autotune_wcc()`) mentioned in the overview; "Where to go next" article table added;
  canned output regenerated to reflect current (M4/M5) labels.

* **Bug fix**: `print.wcc_res` and `print.wcc_surr` no longer display a spurious
  double colon in the Fisher's Z aggregate label. The fix also renames the label from
  `Mean |Fisher's Z|` to `Mean Abs. Fisher's Z` to avoid the cli markup ambiguity.

## M6 — Parameter guidance & synchrony multiverse

### Phase B — autotune_wcc rewrite + vignette

* **`autotune_wcc()` rewritten** as a thin wrapper over `synchrony_multiverse()`
  + a gated stability-penalized selection rule. Takes a `dyad_list`, sweeps a
  seconds-specified grid, and selects the parameter cell that is (a) significant
  vs. the matched-null surrogate in at least `sig_pct` (default 50%) of dyads,
  and (b) maximizes `median(ES) - iqr_penalty * IQR(ES)` across dyads. Interface
  change: `window_sec`/`lag_sec` replace the old `window_multipliers`/
  `lag_multipliers` arguments.

* **`select_specification()`** — new exported helper that implements the gated
  stability-penalized rule on a list of `bsync_multiverse` objects (one per
  dyad). Can be called directly by advanced users.

* **New vignette** `vignettes/choosing-parameters.Rmd` documents the
  three-tool parameter guidance workflow: `suggest_wcc_params()` for a single
  dyad, `synchrony_multiverse()` for visualizing the specification curve, and
  `autotune_wcc()` for multi-dyad datasets. The older `wcc-params` vignette is
  removed; its content (including the window-overlap discussion) is folded into
  the new guide.

* `glance()` on a `bsync_multiverse` now reports both `n_cells` (total
  specifications in the grid) and `n_valid` (cells that produced a computable
  effect size); `pct_significant` is taken over `n_valid`. This disambiguates
  the grid total from the number of cells actually evaluated.

### Phase A — Synchrony multiverse + suggest_wcc_params rework

* **`synchrony_multiverse()`** — new function that sweeps a seconds-specified
  parameter grid across analytic choices (window size, max lag, increment,
  surrogate method, WCC statistic) and evaluates each specification with a
  matched-null surrogate test (Invariant 2). The headline metric is effect size
  vs. the null, not raw synchrony. Supports all three estimators (`"wcc"`,
  `"wdtw"`, `"wgranger"`). Surrogates are generated once per method and reused
  across every cell sharing that method (efficiency seam). Returns a
  `bsync_multiverse` object (`$grid`, `$settings`, `$robustness`).

* **`plot.bsync_multiverse()`** — Simonsohn-style specification curve: top
  panel shows effect sizes sorted by magnitude with significance highlighting;
  bottom panel is a choice dashboard showing which analytic choices each
  specification used. Uses pure ggplot2 + base `grid` package; no new
  dependencies.

* **Tidy interface for `bsync_multiverse`.** `tidy()` returns the full
  specification grid tibble; `glance()` returns a one-row robustness summary
  (n_cells, significance rate, median ES, IQR, sign-consistency); `as_tibble()`
  aliases `tidy()`.

* **`suggest_wcc_params()` reworked.** New signature takes the actual time
  series (`x`, `y`, `sample_rate`) and estimates the dominant behavioral cycle
  from the signal's own PSD via `evaluate_signal_power()` (pass
  `event_duration_sec` to override with a theoretical estimate). Three hard
  constraints are enforced and reported as warnings: the SUSY lag cap
  (`lag_max <= floor(window_size/2)`), a series-length ceiling
  (`window_size <= floor(n/2)`), and a minimum-samples floor.

## M5 — Shared windowed-surface + surrogate framework + tidy interface

* **Unified `$aggregate` slot.** All three estimators (`wcc()`, `wdtw()`,
  `wgranger()`) now return a named numeric `$aggregate` in place of the old
  per-estimator `$fisher_z` / `$mean_distance` scalar. WCC returns
  `c(mean_abs_z = ...)` or `c(peak = ...)` depending on `statistic`; WDTW
  returns `c(mean_distance = ...)`; Granger returns `c(f_xy = ..., f_yx = ...)`.

* **`bsync_surface` superclass.** All three result objects now inherit
  `"bsync_surface"` in addition to their leaf class (`"wcc_res"`, `"wdtw_res"`,
  `"wgranger_res"`). This enables dispatch of shared methods.

* **Shared infrastructure.**
  - `build_surface_grid()`: single source of truth for grid math and the
    `w_max = window_size - 1` boundary (Invariant 4); replaces per-estimator
    copy-paste.
  - `validate_series()` / `validate_window_params()`: shared input validators.
  - `run_surrogate_engine()`: shared surrogate loop used by all three
    `*_surrogate()` wrappers; accepts a prebuilt grid and surrogate matrix and
    an aggregate-only compute function (no `results_df` on the surrogate path,
    Invariant 7); M6-multiverse seam.
  - `build_surface_heatmap()`: shared heatmap scaffold for `plot.wcc_res` and
    `plot.wdtw_res` (axis labels, time_step scaling, zero-lag line, theme).

* **Tidy interface.** `generics::tidy()`, `generics::glance()`, and
  `tibble::as_tibble()` methods for `bsync_surface` objects:
  - `tidy()` returns one row per cell of `results_df`.
  - `glance()` returns a one-row tibble with aggregate(s) + key settings
    (`n_windows`, `window_size`, `lag_max`, etc.).
  - `as_tibble()` is an alias for `tidy()`.
  - Two `glance()` rows can be bound with `dplyr::bind_rows()` for comparison.

* **External oracle validation (AC0 preflight).** Frozen golden values from
  `dtw` (v1.23-3, symmetric1 step pattern, L1 cost) and `lmtest` (v0.9-40,
  full-series Granger F-statistics) are committed in `test-external-oracle.R`
  and checked on every run without introducing live test dependencies.

* **Bug fix — `wdtw_surrogate(fast_method = TRUE)` window alignment.** The fast
  path now evaluates surrogates over exactly the windows of the observed lagged
  surface (edges reserved at both ends), at lag 0. A regression in an interim
  refactor had it use a lag-free grid that shifted windows past the series end,
  over-counting windows and producing spurious out-of-range `NA`s. The slow
  (default) path was unaffected.

## M4 — Selectable WCC aggregate statistic

* **`wcc()` gains a `statistic` argument** (`"mean_abs_z"` | `"peak"`, default
  `"mean_abs_z"`). `"mean_abs_z"` is the SUSY *mean absolute Z* (Tschacher &
  Meier, 2020) — the mean of |Fisher's Z| over all windows and lags — and
  reproduces the pre-M4 behavior exactly. `"peak"` is the rMEA *best-lag*
  convention (Boker et al., 2002): per window, the maximum |Fisher's Z| across
  lags; then mean across windows.

* **`wcc_surrogate()` gains the same `statistic` argument.** The null
  distribution is always built with the same statistic as the observed value
  (Invariant 2). Pass the same value to both functions.

* **Shared internal helper.** Both functions now call `wcc_aggregate()`, a
  single internal helper, guaranteeing the observed and surrogate aggregates
  are numerically identical.

* **Print methods updated.** `print.wcc_res` and `print.wcc_surr` label the
  aggregate line according to the chosen statistic.

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

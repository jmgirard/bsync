# bsync 0.1.0

First public release.

## Synchrony multiverse and parameter guidance

* **`synchrony_multiverse()`** — new function that sweeps a seconds-specified
  parameter grid across analytic choices (window size, max lag, increment,
  surrogate method, WCC statistic) and evaluates each specification with a
  matched-null surrogate test. The headline metric is effect size vs. the null,
  not raw synchrony. Supports all three estimators (`"wcc"`, `"wdtw"`,
  `"wgranger"`). Surrogates are generated once per method and reused across every
  cell sharing that method (efficiency seam). Returns a `bsync_multiverse` object
  (`$grid`, `$settings`, `$robustness`).

* **`plot.bsync_multiverse()`** — Simonsohn-style specification curve: top
  panel shows effect sizes sorted by magnitude with significance highlighting;
  bottom panel is a choice dashboard showing which analytic choices each
  specification used. Uses pure ggplot2 + base `grid` package; no new
  dependencies.

* **Tidy interface for `bsync_multiverse`.** `tidy()` returns the full
  specification grid tibble; `glance()` returns a one-row robustness summary
  (n_cells, significance rate, median ES, IQR, sign-consistency); `as_tibble()`
  aliases `tidy()`.

* **`suggest_wcc_params()`.** Takes the actual time series (`x`, `y`,
  `sample_rate`) and estimates the dominant behavioral cycle from the signal's
  own PSD via `evaluate_signal_power()` (pass `event_duration_sec` to override
  with a theoretical estimate). Three hard constraints are enforced and reported
  as warnings: the SUSY lag cap (`lag_max <= floor(window_size/2)`), a
  series-length ceiling (`window_size <= floor(n/2)`), and a minimum-samples
  floor.

## Parameter auto-tuning

* **`autotune_wcc()`** — a thin wrapper over `synchrony_multiverse()` plus a
  gated stability-penalized selection rule. Takes a `dyad_list`, sweeps a
  seconds-specified grid, and selects the parameter cell that is (a) significant
  vs. the matched-null surrogate in at least `sig_pct` (default 50%) of dyads,
  and (b) maximizes `median(ES) - iqr_penalty * IQR(ES)` across dyads. It returns
  a classed `bsync_autotune` object with a tidy `print()` method that shows only
  the selected parameters and detectability summary; the per-dyad multiverses
  remain available in `$dyad_multiverses`.

* **`select_specification()`** — new exported helper that implements the gated
  stability-penalized rule on a list of `bsync_multiverse` objects (one per
  dyad). Can be called directly by advanced users.

* `glance()` on a `bsync_multiverse` reports both `n_cells` (total
  specifications in the grid) and `n_valid` (cells that produced a computable
  effect size); `pct_significant` is taken over `n_valid`. This disambiguates
  the grid total from the number of cells actually evaluated.

## Shared windowed-surface framework and tidy interface

* **Unified `$aggregate` slot.** All three estimators (`wcc()`, `wdtw()`,
  `wgranger()`) return a named numeric `$aggregate`. WCC returns
  `c(mean_abs_z = ...)` or `c(peak = ...)` depending on `statistic`; WDTW
  returns `c(mean_distance = ...)`; Granger returns `c(f_xy = ..., f_yx = ...)`.

* **`bsync_surface` superclass.** All three result objects inherit
  `"bsync_surface"` in addition to their leaf class (`"wcc_res"`, `"wdtw_res"`,
  `"wgranger_res"`), enabling dispatch of shared methods.

* **Shared infrastructure.**
  - `build_surface_grid()`: single source of truth for grid math and the
    `w_max = window_size - 1` boundary.
  - `validate_series()` / `validate_window_params()`: shared input validators.
  - `run_surrogate_engine()`: shared surrogate loop used by all three
    `*_surrogate()` wrappers; accepts a prebuilt grid and surrogate matrix and
    an aggregate-only compute function (no `results_df` on the surrogate path),
    and is reused by the multiverse engine.
  - `build_surface_heatmap()`: shared heatmap scaffold for `plot.wcc_res` and
    `plot.wdtw_res` (axis labels, time_step scaling, zero-lag line, theme).

* **Tidy interface.** `generics::tidy()`, `generics::glance()`, and
  `tibble::as_tibble()` methods for `bsync_surface` objects:
  - `tidy()` returns one row per cell of `results_df`.
  - `glance()` returns a one-row tibble with aggregate(s) + key settings
    (`n_windows`, `window_size`, `lag_max`, etc.).
  - `as_tibble()` is an alias for `tidy()`.
  - Two `glance()` rows can be bound with `dplyr::bind_rows()` for comparison.

## Selectable WCC aggregate statistic

* **`wcc()` has a `statistic` argument** (`"mean_abs_z"` | `"peak"`, default
  `"mean_abs_z"`). `"mean_abs_z"` is the SUSY *mean absolute Z* (Tschacher &
  Meier, 2020) — the mean of |Fisher's Z| over all windows and lags. `"peak"` is
  the rMEA *best-lag* convention (Boker et al., 2002): per window, the maximum
  |Fisher's Z| across lags, then mean across windows.

* **`wcc_surrogate()` has the same `statistic` argument.** The null
  distribution is always built with the same statistic as the observed value.
  Pass the same value to both functions.

* **Shared internal helper.** Both functions call `wcc_aggregate()`, a single
  internal helper, guaranteeing the observed and surrogate aggregates are
  numerically identical.

* **Print methods.** `print.wcc_res` and `print.wcc_surr` label the aggregate
  line according to the chosen statistic.

## Surrogate testing

* **External oracle validation.** Frozen golden values from `dtw` (v1.23-3,
  symmetric1 step pattern, L1 cost) and `lmtest` (v0.9-40, full-series Granger
  F-statistics) are committed in `test-external-oracle.R` and checked on every
  run without introducing live test dependencies.

* **Bug fix — `wdtw_surrogate(fast_method = TRUE)` window alignment.** The fast
  path evaluates surrogates over exactly the windows of the observed lagged
  surface (edges reserved at both ends), at lag 0. An interim refactor had it use
  a lag-free grid that shifted windows past the series end, over-counting windows
  and producing spurious out-of-range `NA`s. The slow (default) path was
  unaffected.

## Performance

* **WCC core rewritten to an NA-aware prefix-sum algorithm.** `calc_wcc_cpp()`
  preprocesses six masked prefix-sum arrays per lag in O(n) and evaluates each
  window in O(1), replacing the prior O(w\_max) inner loop. Both `na.rm` modes
  are preserved exactly. Speedup on typical configurations: 5–25× depending on
  `window_size` and `lag_max`; see `bench/RESULTS.md` for measured timings.

* **OpenMP removed.** The prefix-sum speedup makes OpenMP unnecessary for the WCC
  core. `SHLIB_OPENMP_CXXFLAGS` has been stripped from `src/Makevars` and
  `src/Makevars.win`; the package is serial-by-default and fully reproducible on
  all platforms.

## Correctness and robustness

* `na.rm` is honored in WCC. `calc_wcc_cpp()` gains an `na_rm` parameter;
  `na.rm = FALSE` in `wcc()` / `wcc_surrogate()` returns `NA` for any window
  containing a missing value instead of silently using pairwise-complete pairs.

* Window-size semantics fixed to exactly `window_size` samples. `w_max` is set to
  `window_size - 1` at the C++ boundary in all three estimators and surrogate
  wrappers; previously each window spanned `window_size + 1` samples (off-by-one).

* Short-series robustness. All three `create_*_df()` builders and surrogate grid
  builders call `cli::cli_abort()` when the series is too short for the chosen
  `window_size` / `lag_max`; index sequences use `seq_len()` throughout.

* `evaluate_signal_power()` no longer auto-prints its plot. The `plot` argument
  has been removed; call `plot()` on the returned `"signal_power_res"` object
  instead, eliminating the `Rplots.pdf` side effect in non-interactive contexts.

* Condition style unified in `R/impute.R` and `R/surrogate_generation.R`. Legacy
  `stop()` / `warning()` calls replaced with `cli::cli_abort()` / `cli::cli_warn()`.

* `leadership_asymmetry()` docs state the centered sliding-window semantics
  explicitly; `min_valid` is validated as a single positive integer.

## Documentation and messaging

* **`vignette("bsync")`** ("Get started") provides a full end-to-end workflow
  walkthrough using `sim_dyad`, a WCC / WDTW / WGC estimator-choice decision
  table, and a reading map into the six deep-dive articles. This is the
  recommended entry point for new users.

* **Structured pkgdown site**: articles are grouped into "Get started", "Core
  estimators", and "Going deeper" sections; the function reference index is
  organized into eight thematic groups (Estimators, Surrogate testing, Optima &
  leadership, Parameter guidance, Preprocessing x2, Tidy interface, Data).

* **`choosing-parameters` vignette** ("Choosing Analysis Parameters") documents
  the three-tool parameter guidance workflow: `suggest_wcc_params()` for a single
  dyad, `synchrony_multiverse()` for visualizing the specification curve, and
  `autotune_wcc()` for multi-dyad datasets. WDTW and Granger `estimator=` support
  is shown.

* **Cross-linking**: `wdtw-workflow`, `wgranger-workflow`, and
  `determine-downsampling` vignettes carry "See also / Next steps" footers
  consistent with `wcc-workflow`.

* **README**: Granger causality is in the headline scope sentence; the
  parameter-guidance tools (`suggest_wcc_params()`, `synchrony_multiverse()`,
  `autotune_wcc()`) appear in the overview; a "Where to go next" article table is
  included.

* **Bug fix**: `print.wcc_res` and `print.wcc_surr` no longer display a spurious
  double colon in the Fisher's Z aggregate label; the label reads
  `Mean Abs. Fisher's Z`.

* **Runnable examples** on all exported functions, all using the bundled
  `sim_dyad` dataset. Compute-heavy examples (WDTW, surrogate testing, the
  multiverse, and autotune) are wrapped in `\donttest{}` and use modest subsets
  or surrogate counts so they run quickly.

* **Plot axis labels** spell out "Lag (tau)" instead of using the Greek letter,
  so the estimator-surface and optima-overlay plots render on graphics devices
  without UTF-8 support.

* **`suggest_wcc_params()` documentation.** The helper derives its timescale from
  the PSD *power cutoff* (the frequency below which 95% of the signal's power
  lies), not the single largest spectral peak — a deliberate, documented choice
  that is robust to the low-frequency / DC power that dominates raw movement
  spectra.

## CRAN readiness and packaging

* **Build artifacts removed from version control.** `src/*.o`, `src/bsync.{so,dll}`,
  `.DS_Store`, and `tests/testthat/Rplots.pdf` were untracked; `.gitignore` now
  carries `*.o`, `*.so`, `*.dll`, `*.dylib`, and `Rplots.pdf` patterns.

* **Documentation complete.** `@return` added to all 18 exported `print()`,
  `plot()`, and `summary()` methods; `R CMD check --as-cran` reports 0 errors /
  0 warnings / 0 notes.

* **DESCRIPTION.** Description field covers all three windowed estimators (WCC,
  WDTW, Granger), surrogate significance testing (circular-shift and
  phase-randomization generators), optima picking, leadership asymmetry index,
  and the full preprocessing pipeline. `Language: en-US` field added.

* **Spell-check clean.** `inst/WORDLIST` created with domain terms, acronyms, and
  proper names; `spelling::spell_check_package()` returns no errors.

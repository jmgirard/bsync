# bsync (development version)

## Windowed phase synchrony

* `wphase()` and `wphase_surrogate()` take a new `statistic` argument.
  `"mean_plv"` (the default) is the mean PLV over all windows and lags, as
  before. `"peak"` takes the largest PLV across lags in each window and
  averages those per-window values, as `wcc(statistic = "peak")` does. The
  surrogate test summarizes every surrogate with the same statistic, and
  `print()`, `glance()`, and `names(aggregate)` show which one was used.
* New article, `vignette("wphase-workflow")`: band-pass filtering before
  phase extraction, `wphase()`, a circular-shift surrogate test, the
  `"peak"` statistic, optima, and plots. It also shows that for rhythmic
  signals the relative phase tracks the lead–lag better than the lag of the
  highest PLV.

## Segment-shuffling surrogates

* New `generate_surrogate_segment()` cuts a series into segments of
  `segment_size` samples, from sample 1 as in the SUSY package (Tschacher &
  Meier, 2020), and builds each surrogate by putting the segments in a new
  order. Orders are drawn uniformly from all orders other than the original,
  without repeats, so a surrogate can leave a segment in place. This keeps
  the add-one p-value exact. The samples after the last whole segment stay in
  place, and the call says how many there are. Missing values move with
  their segments. With 3 segments or fewer, the call warns that no test can
  reach p <= .05.
* A segment test needs every window to lie inside one segment. For `wcc()`,
  `wdtw()`, and `wphase()`, use `segment_size >= window_size + 2 * lag_max`
  and `window_increment = segment_size`. For `wgranger()`, use
  `segment_size = window_size = window_increment`.
* `synchrony_multiverse()` and `autotune_wcc()` accept
  `surrogate_method = "segment"`. A segment cell uses segments of
  `window_size + 2 * lag_max` samples (`window_size` for Granger) and steps
  its windows one segment at a time, so `increment_pct` does not apply to it.
  A segment cell with fewer than 4 segments is skipped with a warning.

## IAAFT surrogates

* New `generate_surrogate_iaaft()` builds surrogates with the iterative
  amplitude-adjusted Fourier transform (IAAFT) of Schreiber & Schmitz (1996).
  Each surrogate matches both the power spectrum and the value distribution
  of the series. With `match = "spectrum"` (the default), the Fourier
  amplitudes match exactly and the values match closely. With
  `match = "values"`, the values match exactly and the spectrum matches
  closely. The spectrum form is the default because, in a WCC size check,
  the values form rejected a true null in 21 of 200 pairs (10.5%) at a
  nominal 5%, and the spectrum form in 9 of the same 200 pairs (4.5%).
  `max_iter` (default 1000) caps the iterations, and the call warns when a
  surrogate reaches it. A series with `NA` is an error.
* `synchrony_multiverse()` and `autotune_wcc()` accept
  `surrogate_method = "iaaft"`. It needs `y` without missing values, and
  the call stops with an error that says so. `autotune_wcc()` reports IAAFT
  non-convergence once for the whole run.

## Both Granger directions in multiverse summaries and plots

* `synchrony_multiverse(estimator = "wgranger")` now also returns
  `$robustness_yx`, the robustness summary of the y -> x direction. It has the
  same fields as `$robustness`, which stays x -> y.
* `print()`, `summary()`, `glance()`, and `plot()` for `bsync_multiverse`
  objects take a `direction` argument: `"xy"` (x -> y, the default) or `"yx"`
  (y -> x). Before, these methods reported x -> y only. For a Granger result,
  `print()` and `summary()` show a "Direction" line, `glance()` has a
  `direction` column, and the plot title names the direction. `"yx"` on a WCC
  or WDTW result is an error. WCC and WDTW results, printed summaries,
  `glance()` columns, and plots do not change.
* A Granger result made before this change has no `$robustness_yx`. With
  `direction = "yx"`, the methods compute it from the grid's `es_yx` and
  `p_yx` columns.
* `summary()` of a `bsync_multiverse` object prints "ES range: none" when no
  cell has a computable ES, in place of `[Inf, -Inf]` and a warning.

## Add-one surrogate p-values and the p <= .05 rule (changes reported results)

* `wcc_surrogate()`, `wdtw_surrogate()`, `wgranger_surrogate()` (both
  directions), and `wphase_surrogate()` now return the add-one p-value
  (b + 1) / (n + 1) (Phipson & Smyth, 2010) in place of b / n. Here n is the
  number of surrogates and b counts the surrogates at least as extreme as the
  observed statistic. The p-value is never 0, and its smallest value is
  1 / (n + 1). Under a null where the observed data are one more draw from
  the surrogate distribution, the test rejects no more often than its nominal
  level. The b / n form gave 0 when b = 0 and rejected too often. A
  single-dyad pseudo-dyad test (`keep_roles = TRUE`) with N dyads gives p of
  at least 1 / N. So p <= .05 needs at least 20 dyads.
* The rule for a significant surrogate p-value is now p <= .05. The four
  surrogate print methods and the `synchrony_multiverse()` summary and print
  use it. `plot.bsync_multiverse()` and the detectability gate of
  `autotune_wcc()` and `select_specification()` also use it. If n + 1 is a
  multiple of 20, an add-one p-value can be exactly .05 (for example, 4 of
  99 surrogates). Such a result is now significant. A test that rejects at
  p <= .05 keeps its nominal size (Phipson & Smyth, 2010). An NA p-value is
  never significant. The `synchrony_multiverse()` summary counts a computable
  cell with an NA p-value as not significant, and `plot.bsync_multiverse()`
  draws it as not significant. `select_specification()` leaves a dyad
  with an NA p-value out of the significance rate.
* The four surrogate print methods show the p-value rounded to 4 digits in
  fixed notation. If that rounding gives 0, they show `< 0.0001`.
* The four surrogate print methods no longer fail on an NA p-value. They say
  that no significance call was made, and `print()` of a
  `wgranger_surrogate()` result names the direction. When `n_surrogates` is
  below 19, they note that no result can reach p <= .05.
* If `n_surrogates` is below 19, `synchrony_multiverse()` now warns, because
  no cell can then reach p <= .05. `autotune_wcc()` gives this warning
  once per call.

## New estimator: windowed phase synchrony

* **`wphase()`** — new windowed phase-synchrony estimator. Instantaneous phase
  is extracted once per series via the analytic signal (`gsignal::hilbert()`);
  each window-by-lag cell reports the phase-locking value (PLV; Lachaux,
  Rodriguez, Martinerie, & Varela, 1999) and the mean relative phase
  (positive = `x` ahead). Returns the shared `bsync_surface` object with
  `print()`, `summary()`, `plot()` (PLV heatmap), and `tidy()`/`glance()`/
  `as_tibble()` support. Input with missing values aborts with guidance
  (impute or trim first) because the FFT-based analytic signal has no coherent
  per-window NA handling.
* **`wphase_surrogate()`** — matched-null surrogate significance test for
  `wphase()`: the mean-PLV aggregate is computed identically on the observed
  and every surrogate series, with the empirical p-value as its upper tail.
* **`pick_optima()`** and **`leadership_asymmetry()`** now accept the phase
  surface (`wphase_res` / `wphase_optima`), searching for peak phase locking.

## Pseudo-dyad surrogates

* **`generate_surrogate_pseudo()`** — new surrogate generator for one dyad in
  a sample of dyads. Its columns are partner series taken from the other
  dyads (pseudo-dyads), cropped start-aligned to the target's length, ready
  for `wcc_surrogate()`, `wdtw_surrogate()`, `wgranger_surrogate()`, and
  `wphase_surrogate()`. This null keeps co-movement that a shared task
  produces and removes coupling specific to the real interaction (the rMEA
  pseudo-synchrony approach; Kleinbub & Ramseyer, 2020). Partners shorter
  than the target are excluded with a warning.
* **`generate_pseudo_dyads()`** — builds the sample-wide set of pseudo-dyads
  (every pairing of series from different dyads, or a random subset), each in
  the `list(x, y)` form, for comparing real dyads to pseudo-dyads. Both
  generators keep partner roles by default (`keep_roles = TRUE`), because
  the lead-lag sign depends on which series is `x`.
* The surrogate-testing vignette has a new pseudo-dyad section with a worked
  example.

## Stricter `dyad_list` reading

* `autotune_wcc()`, `generate_surrogate_pseudo()`, and
  `generate_pseudo_dyads()` now read each dyad by one rule. If its names
  contain `x` and `y` exactly once each, those elements are read by name, in
  any position. If its names contain neither, the first two elements are
  read. Any other use of the names `x` and `y` is an error. Names are no
  longer matched partially, and every error about a dyad's names or shape
  names the dyad's index in `dyad_list`.
* Results change without an error for two kinds of input. A data frame with
  columns named `x` and `y` that are not its first two columns, in that
  order (for example `y, x` or `time, x, y`), is now read by name. Data
  frames used to be read by position. A list whose names only start with
  `x` or `y`, for example `list(yy = , xx = )` or
  `list(xval = , other = , yval = )`, used to be matched partially and is
  now read by position.
* `autotune_wcc()` reads every dyad and checks its names and shape before
  it samples `n_tune_dyads` dyads. A dyad with an unusable naming or shape
  is now an error even if the sample leaves it out. A single data frame
  passed as `dyad_list` is now an error that says to wrap it in a list.
* If `dyad_list` is named, `generate_surrogate_pseudo()` and
  `generate_pseudo_dyads()` carry the names into their output: a `name`
  column (`x_name` and `y_name` for pseudo-dyads) in the `"sources"`
  attribute, plus matrix column names and list element names. The names
  must be non-empty, not `NA`, and unique. Without names, the name columns are `NA`.

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

* **Circular-inference guardrails on `autotune_wcc()`.** Because auto-tuning
  *selects* the parameters that maximize the effect size, re-testing synchrony
  with those parameters on the same dyads is circular (double-dipping). The help
  page now documents this explicitly and points to the two defensible remedies
  (split-sample tuning, or reporting multiverse robustness instead of the tuned
  winner), and `print.bsync_autotune()` emits a matching caution.

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

* **Phase-randomized surrogates preserve the sign of the mean.**
  `generate_surrogate_phase()` previously took the modulus of the DC (and
  even-length Nyquist) Fourier terms, which flipped the sign of the mean for
  negative-mean signals in every surrogate. These real-valued terms are now
  preserved exactly, so the surrogate mean and variance match the original for
  any signal. This is inert for the mean-invariant correlation statistics but
  matters for `wdtw()` with `scale_method = "none"`.

* **Phase-randomized surrogates support odd-length series natively.**
  `generate_surrogate_phase()` no longer aborts on an odd number of
  observations; it constructs a valid conjugate-symmetric spectrum for both
  parities and always returns a matrix with `length(y)` rows. This removes a
  crash in `synchrony_multiverse()` and `autotune_wcc()` (whose default
  `surrogate_method = "phase"` previously errored on any odd-length input). The
  `trim_odd` argument is retained for backward compatibility.

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

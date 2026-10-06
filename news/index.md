# Changelog

## bsync (development version)

### Add-one surrogate p-values (changes reported results)

- [`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md),
  [`wdtw_surrogate()`](https://jmgirard.github.io/bsync/reference/wdtw_surrogate.md),
  [`wgranger_surrogate()`](https://jmgirard.github.io/bsync/reference/wgranger_surrogate.md)
  (both directions), and
  [`wphase_surrogate()`](https://jmgirard.github.io/bsync/reference/wphase_surrogate.md)
  now return the add-one p-value (b + 1) / (n + 1) (Phipson &
  Smyth, 2010) in place of b / n. Here n is the number of surrogates and
  b counts the surrogates at least as extreme as the observed statistic.
  The p-value is never 0, and its smallest value is 1 / (n + 1). Under a
  null where the observed data are one more draw from the surrogate
  distribution, the test rejects no more often than its nominal level.
  The b / n form gave 0 when b = 0 and rejected too often. A single-dyad
  pseudo-dyad test (`keep_roles = TRUE`) with N dyads gives p of at
  least 1 / N. So p \< .05 needs at least 21 dyads.
- The four surrogate print methods show the p-value rounded to 4 digits
  in fixed notation. If that rounding gives 0, they show `< 0.0001`.
- If `n_surrogates < 20`,
  [`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md)
  now warns, because no cell can then reach p \< .05.
  [`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md)
  gives this warning once per call.

### New estimator: windowed phase synchrony

- **[`wphase()`](https://jmgirard.github.io/bsync/reference/wphase.md)**
  — new windowed phase-synchrony estimator. Instantaneous phase is
  extracted once per series via the analytic signal
  ([`gsignal::hilbert()`](https://rdrr.io/pkg/gsignal/man/hilbert.html));
  each window-by-lag cell reports the phase-locking value (PLV; Lachaux,
  Rodriguez, Martinerie, & Varela, 1999) and the mean relative phase
  (positive = `x` ahead). Returns the shared `bsync_surface` object with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) (PLV
  heatmap), and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html)/[`glance()`](https://generics.r-lib.org/reference/glance.html)/
  [`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
  support. Input with missing values aborts with guidance (impute or
  trim first) because the FFT-based analytic signal has no coherent
  per-window NA handling.
- **[`wphase_surrogate()`](https://jmgirard.github.io/bsync/reference/wphase_surrogate.md)**
  — matched-null surrogate significance test for
  [`wphase()`](https://jmgirard.github.io/bsync/reference/wphase.md):
  the mean-PLV aggregate is computed identically on the observed and
  every surrogate series, with the empirical p-value as its upper tail.
- **[`pick_optima()`](https://jmgirard.github.io/bsync/reference/pick_optima.md)**
  and
  **[`leadership_asymmetry()`](https://jmgirard.github.io/bsync/reference/leadership_asymmetry.md)**
  now accept the phase surface (`wphase_res` / `wphase_optima`),
  searching for peak phase locking.

### Pseudo-dyad surrogates

- **[`generate_surrogate_pseudo()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_pseudo.md)**
  — new surrogate generator for one dyad in a sample of dyads. Its
  columns are partner series taken from the other dyads (pseudo-dyads),
  cropped start-aligned to the target’s length, ready for
  [`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md),
  [`wdtw_surrogate()`](https://jmgirard.github.io/bsync/reference/wdtw_surrogate.md),
  [`wgranger_surrogate()`](https://jmgirard.github.io/bsync/reference/wgranger_surrogate.md),
  and
  [`wphase_surrogate()`](https://jmgirard.github.io/bsync/reference/wphase_surrogate.md).
  This null keeps co-movement that a shared task produces and removes
  coupling specific to the real interaction (the rMEA pseudo-synchrony
  approach; Kleinbub & Ramseyer, 2020). Partners shorter than the target
  are excluded with a warning.
- **[`generate_pseudo_dyads()`](https://jmgirard.github.io/bsync/reference/generate_pseudo_dyads.md)**
  — builds the sample-wide set of pseudo-dyads (every pairing of series
  from different dyads, or a random subset), each in the `list(x, y)`
  form, for comparing real dyads to pseudo-dyads. Both generators keep
  partner roles by default (`keep_roles = TRUE`), because the lead-lag
  sign depends on which series is `x`.
- The surrogate-testing vignette has a new pseudo-dyad section with a
  worked example.

### Stricter `dyad_list` reading

- [`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md),
  [`generate_surrogate_pseudo()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_pseudo.md),
  and
  [`generate_pseudo_dyads()`](https://jmgirard.github.io/bsync/reference/generate_pseudo_dyads.md)
  now read each dyad by one rule. If its names contain `x` and `y`
  exactly once each, those elements are read by name, in any position.
  If its names contain neither, the first two elements are read. Any
  other use of the names `x` and `y` is an error. Names are no longer
  matched partially, and every error about a dyad’s names or shape names
  the dyad’s index in `dyad_list`.
- Results change without an error for two kinds of input. A data frame
  with columns named `x` and `y` that are not its first two columns, in
  that order (for example `y, x` or `time, x, y`), is now read by name.
  Data frames used to be read by position. A list whose names only start
  with `x` or `y`, for example `list(yy = , xx = )` or
  `list(xval = , other = , yval = )`, used to be matched partially and
  is now read by position.
- [`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md)
  reads every dyad and checks its names and shape before it samples
  `n_tune_dyads` dyads. A dyad with an unusable naming or shape is now
  an error even if the sample leaves it out. A single data frame passed
  as `dyad_list` is now an error that says to wrap it in a list.
- If `dyad_list` is named,
  [`generate_surrogate_pseudo()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_pseudo.md)
  and
  [`generate_pseudo_dyads()`](https://jmgirard.github.io/bsync/reference/generate_pseudo_dyads.md)
  carry the names into their output: a `name` column (`x_name` and
  `y_name` for pseudo-dyads) in the `"sources"` attribute, plus matrix
  column names and list element names. The names must be non-empty, not
  `NA`, and unique. Without names, the name columns are `NA`.

## bsync 0.1.0

First public release.

### Synchrony multiverse and parameter guidance

- **[`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md)**
  — new function that sweeps a seconds-specified parameter grid across
  analytic choices (window size, max lag, increment, surrogate method,
  WCC statistic) and evaluates each specification with a matched-null
  surrogate test. The headline metric is effect size vs. the null, not
  raw synchrony. Supports all three estimators (`"wcc"`, `"wdtw"`,
  `"wgranger"`). Surrogates are generated once per method and reused
  across every cell sharing that method (efficiency seam). Returns a
  `bsync_multiverse` object (`$grid`, `$settings`, `$robustness`).

- **[`plot.bsync_multiverse()`](https://jmgirard.github.io/bsync/reference/plot.bsync_multiverse.md)**
  — Simonsohn-style specification curve: top panel shows effect sizes
  sorted by magnitude with significance highlighting; bottom panel is a
  choice dashboard showing which analytic choices each specification
  used. Uses pure ggplot2 + base `grid` package; no new dependencies.

- **Tidy interface for `bsync_multiverse`.**
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
  full specification grid tibble;
  [`glance()`](https://generics.r-lib.org/reference/glance.html) returns
  a one-row robustness summary (n_cells, significance rate, median ES,
  IQR, sign-consistency);
  [`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
  aliases [`tidy()`](https://generics.r-lib.org/reference/tidy.html).

- **[`suggest_wcc_params()`](https://jmgirard.github.io/bsync/reference/suggest_wcc_params.md).**
  Takes the actual time series (`x`, `y`, `sample_rate`) and estimates
  the dominant behavioral cycle from the signal’s own PSD via
  [`evaluate_signal_power()`](https://jmgirard.github.io/bsync/reference/evaluate_signal_power.md)
  (pass `event_duration_sec` to override with a theoretical estimate).
  Three hard constraints are enforced and reported as warnings: the SUSY
  lag cap (`lag_max <= floor(window_size/2)`), a series-length ceiling
  (`window_size <= floor(n/2)`), and a minimum-samples floor.

### Parameter auto-tuning

- **[`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md)**
  — a thin wrapper over
  [`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md)
  plus a gated stability-penalized selection rule. Takes a `dyad_list`,
  sweeps a seconds-specified grid, and selects the parameter cell that
  is (a) significant vs. the matched-null surrogate in at least
  `sig_pct` (default 50%) of dyads, and (b) maximizes
  `median(ES) - iqr_penalty * IQR(ES)` across dyads. It returns a
  classed `bsync_autotune` object with a tidy
  [`print()`](https://rdrr.io/r/base/print.html) method that shows only
  the selected parameters and detectability summary; the per-dyad
  multiverses remain available in `$dyad_multiverses`.

- **[`select_specification()`](https://jmgirard.github.io/bsync/reference/select_specification.md)**
  — new exported helper that implements the gated stability-penalized
  rule on a list of `bsync_multiverse` objects (one per dyad). Can be
  called directly by advanced users.

- **Circular-inference guardrails on
  [`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md).**
  Because auto-tuning *selects* the parameters that maximize the effect
  size, re-testing synchrony with those parameters on the same dyads is
  circular (double-dipping). The help page now documents this explicitly
  and points to the two defensible remedies (split-sample tuning, or
  reporting multiverse robustness instead of the tuned winner), and
  [`print.bsync_autotune()`](https://jmgirard.github.io/bsync/reference/print.bsync_autotune.md)
  emits a matching caution.

- [`glance()`](https://generics.r-lib.org/reference/glance.html) on a
  `bsync_multiverse` reports both `n_cells` (total specifications in the
  grid) and `n_valid` (cells that produced a computable effect size);
  `pct_significant` is taken over `n_valid`. This disambiguates the grid
  total from the number of cells actually evaluated.

### Shared windowed-surface framework and tidy interface

- **Unified `$aggregate` slot.** All three estimators
  ([`wcc()`](https://jmgirard.github.io/bsync/reference/wcc.md),
  [`wdtw()`](https://jmgirard.github.io/bsync/reference/wdtw.md),
  [`wgranger()`](https://jmgirard.github.io/bsync/reference/wgranger.md))
  return a named numeric `$aggregate`. WCC returns `c(mean_abs_z = ...)`
  or `c(peak = ...)` depending on `statistic`; WDTW returns
  `c(mean_distance = ...)`; Granger returns `c(f_xy = ..., f_yx = ...)`.

- **`bsync_surface` superclass.** All three result objects inherit
  `"bsync_surface"` in addition to their leaf class (`"wcc_res"`,
  `"wdtw_res"`, `"wgranger_res"`), enabling dispatch of shared methods.

- **Shared infrastructure.**

  - `build_surface_grid()`: single source of truth for grid math and the
    `w_max = window_size - 1` boundary.
  - `validate_series()` / `validate_window_params()`: shared input
    validators.
  - `run_surrogate_engine()`: shared surrogate loop used by all three
    `*_surrogate()` wrappers; accepts a prebuilt grid and surrogate
    matrix and an aggregate-only compute function (no `results_df` on
    the surrogate path), and is reused by the multiverse engine.
  - `build_surface_heatmap()`: shared heatmap scaffold for
    `plot.wcc_res` and `plot.wdtw_res` (axis labels, time_step scaling,
    zero-lag line, theme).

- **Tidy interface.**
  [`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`generics::glance()`](https://generics.r-lib.org/reference/glance.html),
  and
  [`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
  methods for `bsync_surface` objects:

  - [`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns
    one row per cell of `results_df`.
  - [`glance()`](https://generics.r-lib.org/reference/glance.html)
    returns a one-row tibble with aggregate(s) + key settings
    (`n_windows`, `window_size`, `lag_max`, etc.).
  - [`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
    is an alias for
    [`tidy()`](https://generics.r-lib.org/reference/tidy.html).
  - Two [`glance()`](https://generics.r-lib.org/reference/glance.html)
    rows can be bound with
    [`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
    for comparison.

### Selectable WCC aggregate statistic

- **[`wcc()`](https://jmgirard.github.io/bsync/reference/wcc.md) has a
  `statistic` argument** (`"mean_abs_z"` \| `"peak"`, default
  `"mean_abs_z"`). `"mean_abs_z"` is the SUSY *mean absolute Z*
  (Tschacher & Meier, 2020) — the mean of \|Fisher’s Z\| over all
  windows and lags. `"peak"` is the rMEA *best-lag* convention (Boker et
  al., 2002): per window, the maximum \|Fisher’s Z\| across lags, then
  mean across windows.

- **[`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md)
  has the same `statistic` argument.** The null distribution is always
  built with the same statistic as the observed value. Pass the same
  value to both functions.

- **Shared internal helper.** Both functions call `wcc_aggregate()`, a
  single internal helper, guaranteeing the observed and surrogate
  aggregates are numerically identical.

- **Print methods.** `print.wcc_res` and `print.wcc_surr` label the
  aggregate line according to the chosen statistic.

### Surrogate testing

- **External oracle validation.** Frozen golden values from `dtw`
  (v1.23-3, symmetric1 step pattern, L1 cost) and `lmtest` (v0.9-40,
  full-series Granger F-statistics) are committed in
  `test-external-oracle.R` and checked on every run without introducing
  live test dependencies.

- **Bug fix — `wdtw_surrogate(fast_method = TRUE)` window alignment.**
  The fast path evaluates surrogates over exactly the windows of the
  observed lagged surface (edges reserved at both ends), at lag 0. An
  interim refactor had it use a lag-free grid that shifted windows past
  the series end, over-counting windows and producing spurious
  out-of-range `NA`s. The slow (default) path was unaffected.

### Performance

- **WCC core rewritten to an NA-aware prefix-sum algorithm.**
  `calc_wcc_cpp()` preprocesses six masked prefix-sum arrays per lag in
  O(n) and evaluates each window in O(1), replacing the prior O(w_max)
  inner loop. Both `na.rm` modes are preserved exactly. Speedup on
  typical configurations: 5–25× depending on `window_size` and
  `lag_max`; see `bench/RESULTS.md` for measured timings.

- **OpenMP removed.** The prefix-sum speedup makes OpenMP unnecessary
  for the WCC core. `SHLIB_OPENMP_CXXFLAGS` has been stripped from
  `src/Makevars` and `src/Makevars.win`; the package is
  serial-by-default and fully reproducible on all platforms.

### Correctness and robustness

- **Phase-randomized surrogates preserve the sign of the mean.**
  [`generate_surrogate_phase()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_phase.md)
  previously took the modulus of the DC (and even-length Nyquist)
  Fourier terms, which flipped the sign of the mean for negative-mean
  signals in every surrogate. These real-valued terms are now preserved
  exactly, so the surrogate mean and variance match the original for any
  signal. This is inert for the mean-invariant correlation statistics
  but matters for
  [`wdtw()`](https://jmgirard.github.io/bsync/reference/wdtw.md) with
  `scale_method = "none"`.

- **Phase-randomized surrogates support odd-length series natively.**
  [`generate_surrogate_phase()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_phase.md)
  no longer aborts on an odd number of observations; it constructs a
  valid conjugate-symmetric spectrum for both parities and always
  returns a matrix with `length(y)` rows. This removes a crash in
  [`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md)
  and
  [`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md)
  (whose default `surrogate_method = "phase"` previously errored on any
  odd-length input). The `trim_odd` argument is retained for backward
  compatibility.

- `na.rm` is honored in WCC. `calc_wcc_cpp()` gains an `na_rm`
  parameter; `na.rm = FALSE` in
  [`wcc()`](https://jmgirard.github.io/bsync/reference/wcc.md) /
  [`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md)
  returns `NA` for any window containing a missing value instead of
  silently using pairwise-complete pairs.

- Window-size semantics fixed to exactly `window_size` samples. `w_max`
  is set to `window_size - 1` at the C++ boundary in all three
  estimators and surrogate wrappers; previously each window spanned
  `window_size + 1` samples (off-by-one).

- Short-series robustness. All three `create_*_df()` builders and
  surrogate grid builders call
  [`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)
  when the series is too short for the chosen `window_size` / `lag_max`;
  index sequences use [`seq_len()`](https://rdrr.io/r/base/seq.html)
  throughout.

- [`evaluate_signal_power()`](https://jmgirard.github.io/bsync/reference/evaluate_signal_power.md)
  no longer auto-prints its plot. The `plot` argument has been removed;
  call [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the
  returned `"signal_power_res"` object instead, eliminating the
  `Rplots.pdf` side effect in non-interactive contexts.

- Condition style unified in `R/impute.R` and
  `R/surrogate_generation.R`. Legacy
  [`stop()`](https://rdrr.io/r/base/stop.html) /
  [`warning()`](https://rdrr.io/r/base/warning.html) calls replaced with
  [`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html) /
  [`cli::cli_warn()`](https://cli.r-lib.org/reference/cli_abort.html).

- [`leadership_asymmetry()`](https://jmgirard.github.io/bsync/reference/leadership_asymmetry.md)
  docs state the centered sliding-window semantics explicitly;
  `min_valid` is validated as a single positive integer.

### Documentation and messaging

- **[`vignette("bsync")`](https://jmgirard.github.io/bsync/articles/bsync.md)**
  (“Get started”) provides a full end-to-end workflow walkthrough using
  `sim_dyad`, a WCC / WDTW / WGC estimator-choice decision table, and a
  reading map into the six deep-dive articles. This is the recommended
  entry point for new users.

- **Structured pkgdown site**: articles are grouped into “Get started”,
  “Core estimators”, and “Going deeper” sections; the function reference
  index is organized into eight thematic groups (Estimators, Surrogate
  testing, Optima & leadership, Parameter guidance, Preprocessing x2,
  Tidy interface, Data).

- **`choosing-parameters` vignette** (“Choosing Analysis Parameters”)
  documents the three-tool parameter guidance workflow:
  [`suggest_wcc_params()`](https://jmgirard.github.io/bsync/reference/suggest_wcc_params.md)
  for a single dyad,
  [`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md)
  for visualizing the specification curve, and
  [`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md)
  for multi-dyad datasets. WDTW and Granger `estimator=` support is
  shown.

- **Cross-linking**: `wdtw-workflow`, `wgranger-workflow`, and
  `determine-downsampling` vignettes carry “See also / Next steps”
  footers consistent with `wcc-workflow`.

- **README**: Granger causality is in the headline scope sentence; the
  parameter-guidance tools
  ([`suggest_wcc_params()`](https://jmgirard.github.io/bsync/reference/suggest_wcc_params.md),
  [`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md),
  [`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md))
  appear in the overview; a “Where to go next” article table is
  included.

- **Bug fix**: `print.wcc_res` and `print.wcc_surr` no longer display a
  spurious double colon in the Fisher’s Z aggregate label; the label
  reads `Mean Abs. Fisher's Z`.

- **Runnable examples** on all exported functions, all using the bundled
  `sim_dyad` dataset. Compute-heavy examples (WDTW, surrogate testing,
  the multiverse, and autotune) are wrapped in `\donttest{}` and use
  modest subsets or surrogate counts so they run quickly.

- **Plot axis labels** spell out “Lag (tau)” instead of using the Greek
  letter, so the estimator-surface and optima-overlay plots render on
  graphics devices without UTF-8 support.

- **[`suggest_wcc_params()`](https://jmgirard.github.io/bsync/reference/suggest_wcc_params.md)
  documentation.** The helper derives its timescale from the PSD *power
  cutoff* (the frequency below which 95% of the signal’s power lies),
  not the single largest spectral peak — a deliberate, documented choice
  that is robust to the low-frequency / DC power that dominates raw
  movement spectra.

### CRAN readiness and packaging

- **Build artifacts removed from version control.** `src/*.o`,
  `src/bsync.{so,dll}`, `.DS_Store`, and `tests/testthat/Rplots.pdf`
  were untracked; `.gitignore` now carries `*.o`, `*.so`, `*.dll`,
  `*.dylib`, and `Rplots.pdf` patterns.

- **Documentation complete.** `@return` added to all 18 exported
  [`print()`](https://rdrr.io/r/base/print.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods;
  `R CMD check --as-cran` reports 0 errors / 0 warnings / 0 notes.

- **DESCRIPTION.** Description field covers all three windowed
  estimators (WCC, WDTW, Granger), surrogate significance testing
  (circular-shift and phase-randomization generators), optima picking,
  leadership asymmetry index, and the full preprocessing pipeline.
  `Language: en-US` field added.

- **Spell-check clean.** `inst/WORDLIST` created with domain terms,
  acronyms, and proper names;
  [`spelling::spell_check_package()`](https://docs.ropensci.org/spelling//reference/spell_check_package.html)
  returns no errors.

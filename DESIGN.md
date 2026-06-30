# DESIGN.md — `bsync`: Behavioral Synchrony Analyses (R package)

> Working design brief and source of truth for design decisions. It records what the package is,
> the contracts each estimator implements, the recommended defaults and their rationale, and the
> milestone roadmap. `CLAUDE.md` (repo root) covers *how we work*; this file covers *what we're
> building*. When the two disagree, **this file wins for design**, `CLAUDE.md` wins for process —
> flag the conflict.

Package name: **`bsync`**. Main calls: `wcc()`, `wdtw()`, `wgranger()`, plus the preprocessing,
surrogate, optima, and tuning helpers. Site: <https://jmgirard.github.io/bsync/>.

---

## 1. Purpose

Provide a fast, modern, and *educational* toolkit for quantifying **interpersonal / behavioral
synchrony** in continuous dyadic time series (motion capture, computer-vision tracking,
physiology, acoustic features). The package couples efficient, C++-backed estimators of
nonstationary lead–lag structure with (a) a principled preprocessing pipeline, (b) rigorous
null-hypothesis testing via surrogate (pseudo-synchrony) data, and (c) helper functions and
vignettes that teach the user to make good, data-driven choices about downsampling, smoothing, and
hyperparameter selection. It must be **accessible to newcomers** (safe defaults, guidance, clear
output) and **powerful for experts** (every consequential choice is exposed and tunable).

## 2. Scope & positioning

**Prior art (R ecosystem).** The closest siblings, verified, and what they do / don't do:

- **`SUSY`** — *Surrogate Synchrony* (Tschacher & Meier 2020). The canonical implementation of the
  exact method behind bsync's WCC path: windowed cross-correlation of a dyadic (2-column) series,
  Fisher-Z transform of each correlation, **mean Z** and **mean absolute Z** aggregated within and
  then across segments, with an **effect size (ES)** computed against **segment-shuffling
  surrogates** (`susy(x, segment, Hz, maxlag = 3, surrogates.total = 500)`). Pure R; **WCC only** —
  no DTW, no Granger, no peak/valley optima, no preprocessing or tuning layer.
- **`mvSUSY`** — *Multivariate Surrogate Synchrony* (Meier & Tschacher 2021). Extends SUSY beyond
  the dyad to >2 simultaneous series with surrogate controls. The reference for bsync's group /
  multivariate direction (§15, M9).
- **`rMEA`** — *Motion Energy Analysis synchrony* (Kleinbub & Ramseyer 2020). Import/filter/plot MEA
  time series; `MEAccf()` does windowed lagged cross-correlation with increments, a per-window
  **best lag**, pseudosynchrony via shuffled-surrogate dyads, and **lead/follow** indices
  (`s1_lead`/`s2_lead`). Pure R; the MEA-specific sibling to bsync's WCC + optima + leadership layer.
- **`crqa`** (Coco & Dale 2014) — (cross-)recurrence quantification analysis. Reference for M8.
- **`dtw`** (Giorgino 2009) — full-series dynamic time warping (not windowed/nonstationary).
- **`WaveletComp` / `biwavelet`** — wavelet coherence. Reference for M7.
- base `ccf`, `signal`/`gsignal` — cross-correlation and signal-processing primitives.

**Our contribution.** bsync is, in effect, the **SUSY/rMEA windowed-cross-correlation-plus-surrogate
method generalized and accelerated**, plus more:
- **C++/Armadillo cores** for the inner loops, where SUSY/rMEA/crqa are pure R.
- **Multiple estimators behind one surface contract** (§4) — WCC *and* windowed DTW *and* windowed
  Granger — not WCC alone, with a shared optima/surrogate/plot layer.
- **Matched-null surrogate testing** for every estimator, with circular-shift and
  phase-randomization generators (segment-shuffling and IAAFT planned, §6) rather than a single
  shuffle scheme.
- **Peak/valley optima + leadership asymmetry** as a shared, estimator-agnostic layer.
- **An educational preprocessing + tuning layer** — PSD downsampling guidance, zero-phase smoothing,
  kinematics, gap diagnosis/imputation; `suggest_wcc_params()` (theory-based starting points) and
  `autotune_wcc()` (empirical surrogate-driven grid search) — which the siblings do not provide.

**Lineage of the WCC statistics (cite, don't reinvent).** bsync's default WCC aggregate
(`mean_abs_z`, §9) is precisely SUSY's *mean absolute Z*; the `autotune_wcc()` effect size is SUSY's
*ES*; the M4 peak-per-window option is the rMEA / Boker et al. (2002) *best-lag* convention. Docs
should name this lineage explicitly.

**Honesty caveat to state in docs and `print`.** A synchrony surface is a *descriptive* estimate of
moment-to-moment association, not a causal model. Granger "causality" is predictive precedence in
the Wiener–Granger sense, not mechanism. Surrogate p-values test against a specific null (broken
cross-coupling while preserving marginal structure); state which null each method uses.

## 3. Design philosophy & priorities

Owner-stated priorities and the standing stance on each.

- **Efficiency is a first-class feature — and *this is where `bsync` deliberately diverges from a
  pure-R package.*** The heavy inner loops (windowed correlation, DTW dynamic programming, rolling
  Granger regressions, peak-picking) live in **Rcpp/RcppArmadillo** C++ cores. This is the
  package's core value proposition, not an afterthought. Optimize the cores, but **measure before
  and after** (§13) and never let an optimization change results (Invariant 5).
- **Manageable dependencies.** Compiled deps (`Rcpp`, `RcppArmadillo`) and a small Imports set are
  the cost of admission; do not grow `Imports` casually. Plotting and parallel backends are light
  and already present; flag before adding anything heavy.
- **Best practices.** testthat 3e, roxygen2, pkgdown, `vdiffr` plot snapshots, CI (`R CMD check
  --as-cran`, styler/air, lintr), semantic versioning, numerical regression oracles for the C++
  cores.
- **Safe, reproducible, and *loud* defaults.** (See §9.) When the package makes a consequential
  automatic choice (a suggested window, a downsampling factor, an NA policy, a capped lag), it
  **announces it via cli** so a user who never reads an argument still sees what happened.
- **Accessible *and* powerful.** Newcomers get defaults + guidance + readable `print`/`summary`;
  experts get every hyperparameter, the raw surface, and the surrogate machinery. Don't sacrifice
  one audience for the other.

## 4. Estimators & the shared windowed-surface contract

Three estimators today, each a sliding-window analysis over two equal-length series `x`, `y`. The
**design target (M5)** is that all current and future estimators conform to one *windowed-surface
contract* so the optima, surrogate, and plotting layers are written once.

| Estimator | Function | Core | Metric | Lag dimension |
|---|---|---|---|---|
| Windowed cross-correlation | `wcc()` | `calc_wcc_cpp` | Pearson `r` per (window, lag) | yes (±`lag_max`) |
| Windowed dynamic time warping | `wdtw()` | `calc_wdtw_cpp` | DTW distance per (window, lag) | yes (±`lag_max`) |
| Windowed Granger causality | `wgranger()` | `calc_wgranger_cpp` | `F`/`p` for x→y and y→x per window | no (directional, AR-based) |

### 4.1 Surface contract (target form)

Each estimator returns a list-based S3 object (`wcc_res`, `wdtw_res`, `wgranger_res`) with:

```r
structure(list(
  results_df = <data.frame>,   # long: one row per (window i [, lag tau]); metric column(s)
  settings   = <list>,         # window_size, lag_max, increments, method args, has_time
  <aggregate> = <numeric>,     # one-number summary (fisher_z / mean_distance / ...)
), class = c("<estimator>_res", "list"))
```

- **`results_df` columns.** `i` = window position (a sample index, or the real timestamp when
  `time` is supplied). For lagged estimators, `tau` = lag in samples. Then the metric column(s):
  `wcc` (r), `dtw_dist`, or `f_xy/p_xy/f_yx/p_yx`.
- **Grid construction.** From `window_size`, `lag_max`, `window_increment`, `lag_increment`, the R
  wrapper builds the `(i, tau)` grid, calls the C++ core once, and assembles `results_df`. The grid
  builder is shared logic (`create_*_df`), a prime target for unification in M5.
- **C++ core contract.** Cores take `x, y, i_vals, tau_vals, w_max, <method args>` and return the
  metric vector/frame. Cores assume **clean, validated inputs**, do their own **bounds checking**
  (returning `NA` for out-of-range windows), and never message the user directly (Invariant 1).

### 4.2 Window-size semantics (contract — resolved M1)

A window of `window_size` contains **exactly `window_size` samples**. The R wrapper passes
`w_max = window_size - 1` to the C++ core (whose inclusive `0..w_max` loop then spans exactly
`window_size` points). This is the single source of truth for realized window length (Invariant 4).
*(Pre-M1, the cores used `window_size` as `w_max`, silently yielding `window_size + 1` samples —
fixed in M1.)*

### 4.3 Optima (peak/valley picking)

`pick_optima()` consumes a surface and returns the per-window extremum (peak for WCC, valley for
WDTW) and its lag, via `pick_optima_cpp` (local symmetric-search or global). The picked lag time
series feeds `leadership_asymmetry()` (rolling leader–follower index) and the optima-overlay plot.
This layer is estimator-agnostic by construction and should remain so.

## 5. Preprocessing pipeline (educational core)

Helpers are **advisory and non-destructive** (Invariant 3): they transform only what is asked and
never silently change the basis. Each pairs a *diagnostic* with an *action*.

- **Downsampling.** `evaluate_signal_power()` (Welch PSD → cumulative-power cutoff → recommended
  integer downsample factor + clean target Hz) advises; `downsample_signal()` (regular vectors) and
  `aggregate_by_time()` (irregular/timestamped frames) act. The two action functions are documented
  with an explicit "when to use which" contract.
- **Smoothing.** `smooth_signal()` — Savitzky–Golay (default, zero-phase), centered moving average,
  or zero-phase Butterworth — with optional output clamping. `trim_edges()` removes the boundary
  artifacts that symmetric filters mathematically require.
- **Kinematics.** `calc_velocity_1d()` and `calc_speed_{1,2,3}d()` (central/forward/backward finite
  differences with edge filling).
- **Missing data.** `diagnose_ts_gaps()` (gap profile) → `impute_ts_gaps()` (linear/spline, capped
  by `maxgap`, no extrapolation, carries an imputation-summary attribute).

## 6. Surrogate significance testing (the inferential spine)

Generators produce a null in which genuine cross-coupling is destroyed while each series' own
structure is preserved:

- `generate_surrogate_circular()` — circular shifts (preserve autocorrelation; good for behavioral
  data; can respect `lag_max` to guarantee decoupling).
- `generate_surrogate_phase()` — Fourier phase randomization (preserve power spectrum; good for
  continuous/physiological data; Hermitian-symmetric reconstruction; even-length requirement).

Analysis wrappers (`wcc_surrogate()`, `wdtw_surrogate()`, `wgranger_surrogate()`):

1. Compute the **observed aggregate statistic**.
2. Recompute that **same statistic** on each surrogate column (parallelized via `future.apply`,
   `future.seed = TRUE`).
3. Report the empirical p-value as the tail proportion.

**Invariant 2 (matched null).** The statistic used for the observed value and for every surrogate
must be identical; the p-value is its tail. Changing the observed statistic changes the null in
lockstep. (This is why the selectable WCC statistic — §9, M4 — threads through both `wcc()` and
`wcc_surrogate()`.)

**Planned generators (resolved — to build, not to decide).**
- **IAAFT** (iterative amplitude-adjusted Fourier transform; Schreiber & Schmitz 1996) — matches
  both the amplitude distribution *and* the power spectrum, a stricter null than plain phase
  randomization. Committed to the roadmap (§15).
- **Segment shuffling** — the SUSY/rMEA scheme (recombine non-contemporaneous segments). A natural
  third generator that also makes bsync's nulls directly comparable to those packages. Candidate to
  add alongside IAAFT.

## 7. Result objects, console & tidy representations

- **`print`** — compact cli card per estimator (windows, lags, key settings, aggregate). No matrix
  dumps.
- **`summary`** — distributional detail (quantiles of the metric; for Granger, significant-window
  proportions per direction; for optima, lag directionality / leadership breakdown).
- **`plot`** — surface heatmap (WCC diverging `r`; WDTW sequential distance) with optional zero-lag
  reference; Granger renders rolling F or −log10(p) lines per direction; `plot_optima_overlay()`
  draws picked optima on the surface; `plot.bsync_lai()` draws the leadership index.

**Tidy interface (resolved — to add).** broom-style `tidy()` / `glance()` and `as_tibble()` for the
surface and optima objects, so bsync fits the wider tidy ecosystem (`generics` would move to
Imports). `tidy()` returns the long surface (`i`, `tau`, metric); `glance()` returns one row of
run-level meta (estimator, window/lag settings, aggregate, n windows/lags); `as_tibble()` →
`tidy()`. Built as part of the M5 shared-framework refactor so all estimators get it from one place.

## 8. Visualization

ggplot2 throughout (`scales`, `grDevices` for palettes). Surfaces are tiled heatmaps; time-series
views are lines/steps. Plots are snapshot-tested with `vdiffr`. Color choices use colorblind-safe
diverging (RdBu) / sequential (viridis) palettes. Keep layout (data → geometry) separable from
cosmetic arguments so users can restyle freely.

## 9. Defaults (high-stakes — users will not override these)

Principle: **safe, robust, reproducible, self-disclosing.** Every consequential auto-choice is
announced via cli and documented in roxygen with its rationale.

| Decision | Default | Rationale |
|---|---|---|
| `window_size`, `lag_max` | **required** | force a deliberate, theory-informed choice; never silently pick. `suggest_wcc_params()`/`autotune_wcc()` exist to inform it. |
| window length semantics | **exactly `window_size` samples** | contract (§4.2); `w_max = window_size - 1` at the C++ boundary. |
| `window_increment` / `lag_increment` | `1` / `1` | finest-grained surface by default; widen to trade resolution for speed (Boker et al.). |
| `na.rm` (WCC) | `TRUE` (pairwise within window) | tracking data has dropouts; pairwise keeps usable windows. **Must be honored** — `FALSE` returns `NA` for any window containing `NA` (fixed M1; previously a dead argument). |
| `scale_method` (WDTW) | `"global"` | standardize whole series so distances are comparable across windows; `"local"`/`"none"` for special cases. |
| `distance_metric` (WDTW) | `"L2"` | squared Euclidean local cost; `"L1"` for robustness. |
| `ar_order` (Granger) | `1` | parsimonious; expose for users with a model rationale. Window must keep positive residual df. |
| WCC aggregate statistic | **selectable**, default `"mean_abs_z"` | M4: keep mean \|Fisher-z\| over the surface (= SUSY's *mean absolute Z*; Tschacher & Meier 2020); add peak-per-window (the rMEA / Boker best-lag convention). Threads through observed + surrogate (Inv. 2). |
| surrogate method | none default — user picks | `phase` preserves spectrum (continuous data); `circular` preserves autocorrelation (behavioral data). Documented per data type. |
| `n_surrogates` | `100` | enough to explore; **≥ 1000 advised for reporting** (the `print` method warns below 1000). |
| smoothing | Savitzky–Golay, `sg_order = 3`, `window = 5` | zero-phase, preserves peak shape; order 2 for structural trends, > 3 overfits tracking noise. |
| downsample / aggregate method | `"median"` | robust to single-frame tracking glitches; `"mean"` available. |
| `impute_ts_gaps(maxgap)` | `5` | impute only short gaps; longer gaps stay `NA` (no fabricated structure); never extrapolate edges. |
| `evaluate_signal_power(threshold)` | `0.95` | capture 95% of spectral power below the cutoff before recommending a rate. |

### Documentation standard (owner priority)
- Every default above is documented in roxygen with **why**, not just **what** (e.g.
  `suggest_wcc_params()`'s ≈4-cycles-per-window heuristic must be stated and justified).
- An `@details`/vignette explains the surrogate logic and *which null* each method tests.
- Runnable `@examples` for every exported function (use `sim_dyad`); vignettes reproduce an
  end-to-end workflow.
- `@seealso` cross-links estimators ↔ optima ↔ surrogate ↔ tuning helpers.

## 10. Dependencies

| Tier | Packages | Purpose |
|---|---|---|
| LinkingTo | `Rcpp`, `RcppArmadillo` | compiled cores (the point of the package, §3) |
| Imports | `Rcpp`, `cli`, `dplyr`, `rlang`, `utils` | core + console + grid assembly + validation |
| Imports — parallel | `future.apply` | surrogate / autotune parallelism over replicates |
| Imports — viz | `ggplot2`, `scales`, `grDevices` | surfaces, time-series plots, palettes |
| Imports — signal | `gsignal` | Welch PSD, Savitzky–Golay, Butterworth (replaced `signal`) |
| Suggests | `future`, `knitr`, `rmarkdown`, `testthat (>= 3.0.0)`, `vdiffr` | parallel backend, vignettes, tests, plot snapshots |

**Do not grow `Imports` without flagging it.** Prefer base R + the existing stack over new deps.
Compiled work belongs in the C++ cores, not new R dependencies.

## 11. Reproducibility & parallelism

- Stochastic steps (surrogate generation, autotune dyad sampling) depend on the RNG: respect
  `set.seed()`, use `future.seed = TRUE`, and **never reseed internally** (Invariant 6).
- R-level parallelism is opt-in via a `future::plan()`; functions are plan-agnostic.
- **C-level parallelism (OpenMP) — decision pending (M2).** `src/Makevars` currently *links* OpenMP
  but no `#pragma omp` exists. M2 either (a) adopts OpenMP over the independent window/lag index set
  with a thread argument defaulting to serial (for reproducibility and CRAN), documenting the
  nested-parallelism hazard with `future`, or (b) removes the dead flags. Do not ship linked-but-
  unused flags.

## 12. Result-object & C++ boundary invariants

See `CLAUDE.md` "Invariants" for the authoritative, numbered list. In brief: C++ cores are pure and
validated in R; surrogate nulls match the observed statistic; preprocessing is loud and
non-destructive; window length is a contract; optimizations never change results; stochastics are
reproducible; result objects stay light; supplied `time` maps windows to real timestamps.

## 13. Testing & infrastructure

- testthat 3e; `vdiffr` snapshots for every plot method.
- **Numerical regression oracle (Invariant 5).** Any C++ optimization ships with a test asserting
  output equals the prior implementation within tolerance on `sim_dyad` (and a constructed
  NA-laden case). This is the C++ analogue of an algebra-vs-scores cross-check.
- **Benchmark tracking.** Efficiency milestones record before/after timings via a `bench/` script
  (not a test); acceptance criteria cite measured speedups, not "it's faster."
- Short-series and degenerate-parameter guards tested (the `seq_len`/abort path).
- Surrogate sanity: a coupled signal yields a small p-value; an independent pair does not.
- CI: `R CMD check --as-cran` clean (0/0/0 target); styler/air + lintr clean; no tracked build
  artifacts.

## 14. Decisions resolved & remaining

**Resolved (this design round):**
1. Architecture → **unify all estimators behind the §4 windowed-surface contract**; refactor toward
   a shared grid-builder + surrogate engine + optima/plot layer (target milestone M5).
2. Window-size semantics → **exactly `window_size` samples** via `w_max = window_size - 1` (M1).
3. `na.rm` (WCC) → **honored** (`FALSE` ⇒ NA for any NA-containing window); previously dead (M1).
4. WCC aggregate statistic → **selectable**; keep `mean_abs_z`, add peak-per-window; null matches
   observed (M4).
5. Milestone framing → **split** the hardening work into focused milestones (M1–M4) over one bundle.
6. Prior-art positioning → **verified** (§2): SUSY/mvSUSY (Tschacher & Meier 2020; Meier & Tschacher
   2021), rMEA (Kleinbub & Ramseyer 2020), crqa, dtw, WaveletComp/biwavelet. bsync's WCC + surrogate
   + mean\|Z\| path is the SUSY method generalized and C++-accelerated.
7. **tidy / glance / as_tibble** → **will be added** (§7), built in M5 from the shared framework
   (`generics` → Imports at that point).
8. **IAAFT surrogates** → **will be added** (§6, roadmap M10); segment-shuffling surrogates a
   candidate alongside.

9. **OpenMP** → **removed** (M2). The prefix-sum algorithm eliminates the inner w_max loop that
   motivated OpenMP; the serial implementation already achieves 5–25× speedup over the baseline.
   `SHLIB_OPENMP_CXXFLAGS` stripped from `Makevars`/`Makevars.win`; no `#pragma omp` in any source
   file. Serial-by-default and fully reproducible; revisit only if a future estimator has a
   genuinely parallelisable inner loop that prefix sums cannot collapse.

**Remaining / to resolve at the named milestone:**
- **Shared-surface refactor shape** (M5) — how far to merge Granger (directional, no `tau`) into the
  common contract without contorting it; whether to expose a `bsync_surface` superclass.
- A unified **`bsync_ts` preprocessing object** (§5) — future; specify before building.

## 15. Milestone roadmap

**Hardening cycle (current focus — see `CLAUDE.md`):**

1. **M1 — Correctness & robustness.** Honor `na.rm` (plumb through `calc_wcc_cpp`); reconcile
   window-size semantics (`w_max = window_size - 1`, update `n_r` math + docs); short-series guards
   (`seq_len`, explicit abort) across the three `create_*_df` builders and the three surrogate grid
   builders; stop side-effect plotting in `evaluate_signal_power()`; unify condition style to cli in
   `impute.R`/`surrogate_generation.R`; clarify `leadership_asymmetry()` (sliding-window) and
   `suggest_wcc_params()` (4-cycle heuristic) docs + validate `min_valid`. Tests for each.

2. **M2 — Efficiency.** Prefix-sum WCC core (NA-aware cumulative sums so per-(i, τ) work drops to
   the cross-term), numeric-regression oracle vs. reference on `sim_dyad`. Decide OpenMP: adopt over
   the index set (serial default, thread arg, `_OPENMP`-guarded, future-interaction documented) or
   remove the dead flags. `bench/` script with before/after timings.

3. **M3 — CRAN readiness.** Untrack `src/*.o`, `src/*.{so,dll}`, `**/.DS_Store`,
   `tests/testthat/Rplots.pdf`; extend `.gitignore`/`.Rbuildignore`; regenerate `RcppExports`
   cleanly; `R CMD check --as-cran` to 0/0/0 (examples runnable, no filesystem side effects,
   `\value` everywhere, vignettes build, conditional `Suggests` use); README/pkgdown pass.

4. **M4 — Selectable WCC aggregate statistic.** `statistic = c("mean_abs_z", "peak")` on `wcc()` and
   `wcc_surrogate()`; refactor `fisher_z()` into a dispatch; null matches observed; document in the
   WCC vignette.

**Future methods (later cycles — build on the M5 shared framework):**

5. **M5 — Shared windowed-surface + surrogate framework + tidy interface.** Factor the common
   grid-builder, surrogate engine, optima, and plot layer so new estimators plug in cheaply; resolve
   the Granger-into-the-contract question (§14). Add broom-style `tidy()`/`glance()`/`as_tibble()`
   (§7) for all estimators from this shared layer (`generics` → Imports).
6. **M6 — Phase synchrony (Hilbert).** Analytic-signal instantaneous phase; windowed phase-locking
   value / phase synchrony; relative-phase output.
7. **M7 — Wavelet coherence.** Cross-wavelet / wavelet coherence for nonstationary time–frequency
   lead–lag across scales.
8. **M8 — CRQA / MEA conventions.** Cross-recurrence quantification analysis and motion-energy-
   analysis-style windowed cross-correlation conventions/conveniences.
9. **M9 — Group-level / multivariate workflow.** Tidy multi-dyad pipeline (list / nested frame in,
   per-dyad surfaces out) + aggregation / mixed-model summaries across many dyads, reusing the
   `autotune_wcc()` dyad-list conventions; the `mvSUSY` multivariate-synchrony measures are the
   reference for the >2-series case.
10. **M10 — Expanded surrogate generators.** IAAFT (Schreiber & Schmitz 1996) and segment-shuffling
    (SUSY/rMEA) generators, wired through the matched-null engine (Inv. 2); makes bsync's nulls
    comparable to SUSY/rMEA.

Later, unscheduled: a unified `bsync_ts` preprocessing object, expanded educational vignettes
(choosing a method, interpreting a surface, reporting synchrony).

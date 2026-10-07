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
  multivariate direction (§15, M11).
- **`rMEA`** — *Motion Energy Analysis synchrony* (Kleinbub & Ramseyer 2020). Import/filter/plot MEA
  time series; `MEAccf()` does windowed lagged cross-correlation with increments, a per-window
  **best lag**, pseudosynchrony via shuffled-surrogate dyads, and **lead/follow** indices
  (`s1_lead`/`s2_lead`). Pure R; the MEA-specific sibling to bsync's WCC + optima + leadership layer.
- **`crqa`** (Coco & Dale 2014) — (cross-)recurrence quantification analysis. Reference for M10.
- **`dtw`** (Giorgino 2009) — full-series dynamic time warping (not windowed/nonstationary).
- **`WaveletComp` / `biwavelet`** — wavelet coherence. Reference for M9.
- base `ccf`, `signal`/`gsignal` — cross-correlation and signal-processing primitives.

**Our contribution.** bsync is, in effect, the **SUSY/rMEA windowed-cross-correlation-plus-surrogate
method generalized and accelerated**, plus more:
- **C++/Armadillo cores** for the inner loops, where SUSY/rMEA/crqa are pure R.
- **Multiple estimators behind one surface contract** (§4) — WCC *and* windowed DTW *and* windowed
  Granger — not WCC alone, with a shared optima/surrogate/plot layer.
- **Matched-null surrogate testing** for every estimator, with circular-shift,
  phase-randomization, IAAFT, segment-shuffling (§6; corrected M015), and pseudo-dyad generators rather
  than a single shuffle scheme.
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

Four estimators today (wphase added post-1.0), each a sliding-window analysis over two equal-length series `x`, `y`. The
**design target (M5)** is that all current and future estimators conform to one *windowed-surface
contract* so the optima, surrogate, and plotting layers are written once.

| Estimator | Function | Core | Metric | Lag dimension |
|---|---|---|---|---|
| Windowed cross-correlation | `wcc()` | `calc_wcc_cpp` | Pearson `r` per (window, lag) | yes (±`lag_max`) |
| Windowed dynamic time warping | `wdtw()` | `calc_wdtw_cpp` | DTW distance per (window, lag) | yes (±`lag_max`) |
| Windowed Granger causality | `wgranger()` | `calc_wgranger_cpp` | `F`/`p` for x→y and y→x per window | no (directional, AR-based) |
| Windowed phase synchrony | `wphase()` | `calc_wphase_cpp` | PLV + mean relative phase per (window, lag) | yes (±`lag_max`) |

### 4.1 Surface contract (target form)

Each estimator returns a list-based S3 object (`wcc_res`, `wdtw_res`, `wgranger_res`, `wphase_res`) with:

```r
structure(list(
  results_df = <data.frame>,   # long: one row per (window i [, lag tau]); metric column(s)
  settings   = <list>,         # window_size, lag_max, increments, method args, has_time
  <aggregate> = <numeric>,     # one-number summary (fisher_z / mean_distance / ...)
), class = c("<estimator>_res", "list"))
```

- **`results_df` columns.** `i` = window position (a sample index, or the real timestamp when
  `time` is supplied). For lagged estimators, `tau` = lag in samples. Then the metric column(s):
  `wcc` (r), `dtw_dist`, `f_xy/p_xy/f_yx/p_yx`, or `plv/rel_phase`.
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
  with an explicit "when to use which" contract. *Principle:* choose the target rate from each
  signal's **own informative bandwidth** (the PSD cutoff), at ≥ 2–3× it — not by copying another
  study's rate. Different modalities of the "same" behavior legitimately want different rates (e.g.
  video-derived facial AUs are band-limited well below surface-EMG envelopes), and the PSD makes that
  visible rather than a judgment call. *Planned enhancement (small):* `downsample_signal()` currently
  decimates by **median/mean binning** (robust to tracking glitches but only a partial anti-alias);
  add an **optional Butterworth anti-alias low-pass** before decimation so high-frequency noise is
  not folded into band for non-glitchy signals (default stays median binning; opt-in low-pass).
- **Smoothing.** `smooth_signal()` — Savitzky–Golay (default, zero-phase), centered moving average,
  or zero-phase Butterworth — with optional output clamping. `trim_edges()` removes the boundary
  artifacts that symmetric filters mathematically require.
- **Kinematics.** `calc_velocity_1d()` and `calc_speed_{1,2,3}d()` (central/forward/backward finite
  differences with edge filling).
- **Missing data.** `diagnose_ts_gaps()` (gap profile) → `impute_ts_gaps()` (linear/spline, capped
  by `maxgap`, no extrapolation, carries an imputation-summary attribute).

## 6. Surrogate significance testing (the inferential spine)

Generators produce a null in which genuine cross-coupling is destroyed. The within-dyad
generators preserve each series' own structure. The pseudo-dyad generators replace the partner
with a real series from another dyad (corrected M009):

- `generate_surrogate_circular()` — circular shifts (preserve autocorrelation; good for behavioral
  data; can respect `lag_max` to guarantee decoupling).
- `generate_surrogate_phase()` — Fourier phase randomization (preserve power spectrum; good for
  continuous/physiological data; Hermitian-symmetric reconstruction; even-length requirement).
- `generate_surrogate_iaaft()` — IAAFT (M014; Schreiber & Schmitz 1996,
  `cairn/references/schreiber1996.md`): iterate between the data's Fourier amplitudes and its
  exact values until the rank-ordering stops changing. Each `match` form keeps one property
  exactly and the other closely: `"spectrum"` (default) exact amplitudes, close values;
  `"values"` exact values, close spectrum. The default is the spectrum form because the
  rank-ordered form lowers autocorrelation and made a WCC test liberal (M014 Decisions). `NA`
  aborts.
- `generate_surrogate_segment()` — segment shuffling (M015; SUSY segment cut from sample 1,
  `cairn/references/tschacher2020.md`): each surrogate reorders the k whole segments of `y`
  in an order drawn uniformly from the k! - 1 non-identity orders, without repeats; the tail
  stays in place; `NA` moves with its segment. Exact only when every lagged window lies inside
  one segment: `segment_size >= window_size + 2 * lag_max` with `window_increment =
  segment_size` (wcc/wdtw/wphase), `segment_size = window_size = window_increment` (wgranger);
  wphase stays approximate (global Hilbert transform). The multiverse forces this design for
  `"segment"` cells and skips cells with fewer than 4 segments. D-004; RR01 (archived).
- `generate_surrogate_pseudo()` / `generate_pseudo_dyads()` — pseudo-dyads (M009): partners drawn
  from *other* dyads in a `dyad_list`, cropped start-aligned. A different null from the two above:
  it keeps task-locked co-movement and removes interaction-specific coupling. Convention and
  departures (roles kept by default, no NA-row removal): `cairn/references/kleinbub2020.md`.

Analysis wrappers (`wcc_surrogate()`, `wdtw_surrogate()`, `wgranger_surrogate()`,
`wphase_surrogate()`; corrected M009):

1. Compute the **observed aggregate statistic**.
2. Recompute that **same statistic** on each surrogate column (parallelized via `future.apply`,
   `future.seed = TRUE`).
3. Report the add-one p-value p = (b + 1) / (n + 1). Here n is the number of surrogates, and b
   counts the surrogate statistics at least as extreme as the observed one. The tail is upper for
   wcc, wgranger, and wphase, and lower for wdtw. This is the exact Monte Carlo p-value,
   `cairn/references/phipson2010.md` (p. 6). It is never 0, and its smallest value is
   1 / (n + 1). The b / n form it replaced can be 0 and exceeds the nominal size (p. 4). D-001,
   corrected M011.

**Invariant 2 (matched null).** The statistic used for the observed value and for every surrogate
must be identical; the p-value is its tail, counted in the add-one form of step 3 (phipson2010,
p. 6; corrected M011). Changing the observed statistic changes the null in
lockstep. (This is why the selectable WCC statistic — §9, M4 — threads through both `wcc()` and
`wcc_surrogate()`.)

**Planned generators (resolved — to build, not to decide).** IAAFT shipped in M014 and segment
shuffling in M015 (both listed above; corrected M015). Segment rotations (every segment moved,
exact) are a ROADMAP candidate.

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
| `n_surrogates` | `100` | enough to explore; **≥ 1000 advised for reporting** (the `print` method warns below 1000). `synchrony_multiverse()` and `autotune_wcc()` warn below 19, where the add-one p-value cannot reach p <= .05, the significance rule (D-002), and the surrogate `print` methods note it (corrected M012). |
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
- **Oracle records:** provenance headers in the asserting test files
  (`tests/testthat/test-external-oracle.R`, per-estimator test files) plus the
  committed generators under `data-raw/` (`wphase_null_bound.R`,
  `wphase_mne_pin.py`); each names its source, versions, and asserting test. Surrogate
  generator oracles (M014, M015): `test-surrogate-iaaft.R` and `test-surrogate-segment.R`
  (closed-form references, invariants; the segment file also holds uniformity tests),
  `test-surrogate-calibration.R` (layer-3 size checks with a planted-defect control for
  segments), and `test-multiverse-segment.R` (multiverse cells equal the wrapper calls).
- **Layered validation (four modes, distinct failure types).** Unit tests catch *implementation*
  bugs but cannot catch a *definitional* one: a pure-R oracle written from the same mental model as
  the C++ core agrees with it and is wrong in the same way. The four layers close different gaps:
  1. **Atomic numerical oracle (internal).** Pure-R recomputation of one window/lag against base
     `stats::cor`/`ccf`, matched to ~1e-9 (e.g. the WCC oracle). Catches implementation bugs.
  2. **Convention oracle (external, one-time, *frozen*).** Each *named* method is pinned once,
     during development, against its reference package on a tiny fixed input — WCC/mean\|Z\| vs
     **SUSY**, best-lag vs **rMEA** (`MEAccf`), DTW distance vs **`dtw`** (Giorgino), Granger F vs
     **`lmtest::grangertest`**; future estimators vs their §2 sibling (phase-locking; wavelet
     coherence vs **`biwavelet`/`WaveletComp`**; CRQA vs **`crqa`**). The numbers the external
     package produced are committed as a **static golden fixture** with a provenance comment
     (package + version + exact call); tests run against the frozen vector, **not** the live
     package. **The §2 sibling packages are validation oracles, not just positioning.** Discipline:
     do **not** add these packages as test dependencies (CI/version fragility, convention drift);
     where conventions differ (SUSY segments vs sliding windows), validate the shared kernel
     configured to match. This is the only layer that catches a misunderstanding of the field's
     definition.
  3. **Statistical calibration (simulation, not matching).** The surrogate spine is validated by
     *properties*, not by matching another RNG: under an independent pair the p-value is ~uniform
     (Type-I), under a coupled pair it has power, and each generator preserves its invariant
     (circular shift → autocorrelation; phase randomization → power spectrum, within tolerance).
     Tolerance-banded Monte-Carlo, `skip_on_cran` if slow.
  4. **Regression oracle (Invariant 5).** Any C++ optimization ships with a test asserting output
     equals the prior implementation within tolerance on `sim_dyad` (and a constructed NA-laden
     case). The C++ analogue of an algebra-vs-scores cross-check. **A Layer-2 check must precede a
     Layer-4 freeze** — otherwise a characterization/regression baseline silently locks in a
     definitional bug.
- **Benchmark tracking.** Efficiency milestones record before/after timings via a `bench/` script
  (not a test); acceptance criteria cite measured speedups, not "it's faster."
- Short-series and degenerate-parameter guards tested (the `seq_len`/abort path).
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
8. **IAAFT surrogates** → **added** in M014 (§6; corrected M014), with a `match` argument whose
   default keeps the spectrum exact. Segment-shuffling surrogates are **added** in M015 (§6;
   corrected M015), with uniform non-identity orders and aligned windows (D-004).
   Pseudo-dyad generators shipped in M009.

9. **OpenMP** → **removed** (M2). The prefix-sum algorithm eliminates the inner w_max loop that
   motivated OpenMP; the serial implementation already achieves 5–25× speedup over the baseline.
   `SHLIB_OPENMP_CXXFLAGS` stripped from `Makevars`/`Makevars.win`; no `#pragma omp` in any source
   file. Serial-by-default and fully reproducible; revisit only if a future estimator has a
   genuinely parallelisable inner loop that prefix sums cannot collapse.

10. **Parameter selection** → **resolved (M6).** No single "correct" WCC parameter set exists;
    the optimum is a function of the signal's own timescales. Three read-outs on the shared M5
    engine: (a) `suggest_wcc_params(x, y, sample_rate, ...)` is the PSD-data-driven starting point
    for a single dyad (4-cycles heuristic + three hard SUSY constraints enforced and warned); (b)
    `synchrony_multiverse()` sweeps a seconds-specified grid with matched-null surrogates per cell
    and returns a `bsync_multiverse` object with a Simonsohn-style specification-curve plot; (c)
    `autotune_wcc(dyad_list, ...)` is a thin wrapper over `synchrony_multiverse()` + the gated
    stability-penalized `select_specification()` rule (detectability gate: significant in >= 50%
    of dyads; score: `median(ES) - 0.5 * IQR(ES)`). All three validated against `sim_dyad`;
    `autotune_wcc()` and `suggest_wcc_params()` confirmed on the 0.5 Hz synchrony structure.
    **Efficiency seam (decided M5, exploited M6):** surrogate *generation* stays decoupled from
    surrogate *analysis* (the pre-generated `y_surrogates` matrix interface), and the M5
    `run_surrogate_engine()` accepts a prebuilt grid + prebuilt surrogate matrix and exposes an
    aggregate-only path (core -> aggregate, no per-cell `results_df`). This is what keeps the
    multiverse affordable: surrogates of `y` depend only on `surrogate_method`, not on
    `window_size`/`lag_max`/increment, so one surrogate matrix is generated per method and reused
    across every parameter cell that shares it -- collapsing the dominant cost by the size of the
    parameter grid. The heavy inner compute remains the M2 prefix-sum C++ core; M5/M6 add only R
    orchestration over it, parallelized via `future.apply`.

11. **Shared-surface refactor shape** → **resolved (M5).** Granger joins the *shared machinery*
    (validator, `build_surface_grid()`, `run_surrogate_engine()`, `bsync_surface` superclass, tidy
    layer) but keeps its *estimator-specific surface* — lag-free `results_df` with directional
    `f_xy/p_xy/f_yx/p_yx` and its own line plot — rather than being contorted into a symmetric,
    `tau`-indexed similarity. The **`bsync_surface` superclass** is exposed: `wcc_res`/`wdtw_res`/
    `wgranger_res` inherit it, shared/tidy methods dispatch on it, estimator-specifics on the leaf.
    The per-estimator `$fisher_z`/`$mean_distance` slots were dropped for a uniform named-numeric
    `$aggregate` (`mean_abs_z`/`peak`, `mean_distance`, `c(f_xy, f_yx)`).

**Remaining / to resolve at the named milestone:**
- A unified **`bsync_ts` preprocessing object** (§5) — future; specify before building.

## 15. Milestone roadmap

**The forward roadmap now lives in `cairn/ROADMAP.md`** (the only authority
on milestone status; the M8–M12 estimator ideas are its candidate rows).
Completed milestones M1–M7 have their detailed log entombed verbatim in
`cairn/legacy/MILESTONES.md`; do not duplicate either here. (Replaced at
cairn-init migration, 2026-08-29; the pre-migration §15 text is in git
history.)

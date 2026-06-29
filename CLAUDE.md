# CLAUDE.md — bsync

Operating manual for AI-assisted development of this package. Read
`DESIGN.md` (repo root) first and treat it as the **source of truth**
for all design decisions; this file covers *how we work*, not *what
we’re building*. When this file and `DESIGN.md` disagree, `DESIGN.md`
wins for design and this file wins for process — and flag the conflict.

## What this is

`bsync` is an R package for analyzing **interpersonal / behavioral
synchrony** in continuous dyadic time series. It provides
C++/Armadillo-accelerated windowed estimators of nonstationary lead–lag
structure (windowed cross-correlation, windowed dynamic time warping,
windowed Granger causality), a preprocessing pipeline (PSD-based
downsampling guidance, zero-phase smoothing, kinematics, gap imputation,
time-bin aggregation), surrogate (pseudo-synchrony) significance
testing, peak/valley optima picking and leadership-asymmetry indices,
and theory- and data-driven hyperparameter helpers. Full rationale,
contracts, the surface object spec, and resolved defaults are in
`DESIGN.md`.

The package must serve **newcomers** (safe, loud defaults + guidance)
and **experts** (every hyperparameter exposed, raw surface accessible,
full surrogate machinery) at once.

## Baseline (pre-milestone)

Substantial functionality predates formal milestone tracking; new
milestones begin at **M1**. The baseline, as built:

- **Estimators:**
  [`wcc()`](https://jmgirard.github.io/bsync/reference/wcc.md),
  [`wdtw()`](https://jmgirard.github.io/bsync/reference/wdtw.md),
  [`wgranger()`](https://jmgirard.github.io/bsync/reference/wgranger.md)
  with Rcpp/Armadillo cores (`calc_wcc_cpp`, `calc_wdtw_cpp`,
  `calc_wgranger_cpp`) and S3 `print`/`summary`/`plot`.
- **Optima & leadership:**
  [`pick_optima()`](https://jmgirard.github.io/bsync/reference/pick_optima.md)
  (+ `pick_optima_cpp`),
  [`leadership_asymmetry()`](https://jmgirard.github.io/bsync/reference/leadership_asymmetry.md),
  [`plot_optima_overlay()`](https://jmgirard.github.io/bsync/reference/plot_optima_overlay.md).
- **Surrogates:**
  [`generate_surrogate_circular()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_circular.md),
  [`generate_surrogate_phase()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_phase.md),
  [`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md),
  [`wdtw_surrogate()`](https://jmgirard.github.io/bsync/reference/wdtw_surrogate.md),
  [`wgranger_surrogate()`](https://jmgirard.github.io/bsync/reference/wgranger_surrogate.md)
  (parallel via `future.apply`).
- **Preprocessing:**
  [`evaluate_signal_power()`](https://jmgirard.github.io/bsync/reference/evaluate_signal_power.md),
  [`downsample_signal()`](https://jmgirard.github.io/bsync/reference/downsample_signal.md),
  [`aggregate_by_time()`](https://jmgirard.github.io/bsync/reference/aggregate_by_time.md),
  [`smooth_signal()`](https://jmgirard.github.io/bsync/reference/smooth_signal.md),
  [`trim_edges()`](https://jmgirard.github.io/bsync/reference/trim_edges.md),
  [`calc_velocity_1d()`](https://jmgirard.github.io/bsync/reference/calc_velocity_1d.md),
  `calc_speed_{1,2,3}d()`,
  [`diagnose_ts_gaps()`](https://jmgirard.github.io/bsync/reference/diagnose_ts_gaps.md),
  [`impute_ts_gaps()`](https://jmgirard.github.io/bsync/reference/impute_ts_gaps.md).
- **Tuning:**
  [`suggest_wcc_params()`](https://jmgirard.github.io/bsync/reference/suggest_wcc_params.md),
  [`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md).
- **Data:** `sim_dyad` (simulated 3D dyadic positions with a shifting
  lead–lag relationship).
- **Infra:** testthat 3e suite, `vdiffr` plot snapshots, pkgdown config,
  vignettes (WCC workflow, WDTW, Granger, surrogate testing,
  downsampling, WCC params).

A read of the baseline surfaced the defects M1–M3 address (see Current
focus and `DESIGN.md` §14/§15).

## Completed milestones

*(none yet — M1 is the first formal milestone; see Current focus)*

## Current focus

**Hardening cycle toward a near-term CRAN submission**, run as four
focused milestones via the plan → implement → review loop. **M1 is
active.**

- **M1 — Correctness & robustness (active).** Acceptance criteria:
  1.  `na.rm` is honored in WCC: `calc_wcc_cpp` gains an NA-policy
      parameter; `na.rm = FALSE` returns `NA` for any window containing
      `NA`; forwarded from `create_wcc_df()`. Test proves both modes.
  2.  Window-size semantics fixed to **exactly `window_size` samples**
      via `w_max = window_size - 1` at the C++ boundary in all three
      estimators; `n_r` grid math and all `window_size` docs updated; a
      test asserts realized window length.
  3.  Short-series robustness: `1:n_r` / `0:(n_r-1)` →
      [`seq_len()`](https://rdrr.io/r/base/seq.html); explicit
      [`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)
      when a series is too short for the chosen `window_size`/`lag_max`,
      across the three `create_*_df` builders and the three surrogate
      grid builders. Test for the abort.
  4.  [`evaluate_signal_power()`](https://jmgirard.github.io/bsync/reference/evaluate_signal_power.md)
      no longer auto-[`print()`](https://rdrr.io/r/base/print.html)s its
      plot (returns it); no `Rplots.pdf` side effect; other compute
      functions audited for stray devices.
  5.  Condition style unified to `cli` in `R/impute.R` and
      `R/surrogate_generation.R`.
  6.  Docs:
      [`leadership_asymmetry()`](https://jmgirard.github.io/bsync/reference/leadership_asymmetry.md)
      clarified as a centered **sliding** window + `min_valid`
      validated;
      [`suggest_wcc_params()`](https://jmgirard.github.io/bsync/reference/suggest_wcc_params.md)
      4-cycles-per-window heuristic stated and justified.
- **M2 — Efficiency.** Prefix-sum WCC core (NA-aware) +
  numerical-regression oracle vs. reference on `sim_dyad`; OpenMP
  adopt-or-remove decision (serial default, `_OPENMP`-guarded thread arg
  if adopted); `bench/` before/after timings.
- **M3 — CRAN readiness.** Untrack build artifacts; clean
  `.gitignore`/`.Rbuildignore`; regenerate `RcppExports`;
  `R CMD check --as-cran` → 0/0/0; README/pkgdown pass.
- **M4 — Selectable WCC aggregate statistic.**
  `statistic = c("mean_abs_z","peak")` on
  [`wcc()`](https://jmgirard.github.io/bsync/reference/wcc.md) +
  [`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md);
  null matches observed; vignette documents the choice.

See `DESIGN.md` §15 for the full roadmap (M5 shared framework; M6 phase
synchrony; M7 wavelet coherence; M8 CRQA/MEA; M9 group-level workflow).

## Invariants — do not violate without flagging

These encode hard-won reasoning. Changing them is a design decision, not
a refactor.

1.  **C++ cores are pure and validated in R.** All argument checking, NA
    policy, and grid construction live in the R wrapper; cores assume
    clean inputs, do their own bounds checks (returning `NA` out of
    range), and never message the user. Never expose a `*_cpp` function
    as user API.
2.  **Surrogate nulls match the observed statistic.** The aggregate
    statistic computed on the observed data and on every surrogate must
    be identical; the p-value is its tail. Change one, change both in
    lockstep.
3.  **Loud, non-destructive preprocessing.** Diagnostics advise; action
    functions transform only what is asked and never silently change the
    basis. Announce consequential auto-choices via cli.
4.  **`window_size` means exactly `window_size` samples.** Realized
    window length is a contract; `w_max = window_size - 1` at the C++
    boundary is its single source of truth.
5.  **Optimizations never change results.** Any change to a C++ core
    reproduces the prior implementation’s output within tolerance on
    `sim_dyad` (and an NA-laden case). Speed is free; numbers are
    sacred.
6.  **Reproducible stochastics.** Surrogate generation and autotune
    sampling respect [`set.seed()`](https://rdrr.io/r/base/Random.html)
    and `future.seed = TRUE`; never reseed internally.
7.  **Light result objects.** Carry `results_df` + `settings` + the
    aggregate; raw surrogate draws live only in surrogate objects; raw
    input data is not stored.
8.  **Time integrity.** When `time` is supplied, window positions map to
    real timestamps; edge trimming preserves the true timeline (the
    documented reason `time` exists).

## Resolved defaults (see `DESIGN.md` §9)

`window_size`/`lag_max` **required** · window length = exactly
`window_size` samples · increments = `1`/`1` · WCC `na.rm = TRUE`
(pairwise; honored — `FALSE` ⇒ NA window) · WDTW
`scale_method = "global"`, `distance_metric = "L2"` · Granger
`ar_order = 1` · WCC aggregate statistic selectable, default
`"mean_abs_z"` (M4) · surrogate method user-chosen (phase = spectrum,
circular = autocorrelation) · `n_surrogates = 100` (≥ 1000 advised for
reporting) · smoothing = Savitzky–Golay order 3 · downsample/aggregate =
median · `impute maxgap = 5`, no extrapolation · PSD `threshold = 0.95`.
Don’t change these silently.

## Dependencies (see `DESIGN.md` §10)

Compiled cores via **`LinkingTo: Rcpp, RcppArmadillo`** are the
package’s reason for existing — unlike a pure-R package, heavy inner
loops belong in C++ here. Current `Imports`: `Rcpp`, `cli`, `dplyr`,
`future.apply`, `ggplot2`, `grDevices`, `gsignal`, `rlang`, `scales`,
`utils`. `Suggests`: `future`, `knitr`, `rmarkdown`,
`testthat (>= 3.0.0)`, `vdiffr`. **Do not grow `Imports` without
flagging it** — prefer base R or the existing stack, and put heavy
compute in the C++ cores, not new R deps.

## Dev workflow

R ≥ 4.1 (native pipe `|>` and `\(x)` lambdas allowed). Standard devtools
loop, **with the compiled step**:

``` r

devtools::load_all()              # compiles changed C++ + loads; run after any src/ edit
Rcpp::compileAttributes()         # regenerate RcppExports after changing a [[Rcpp::export]] signature
devtools::document()              # regenerate roxygen docs + NAMESPACE after any roxygen change
devtools::test()                  # run testthat suite
devtools::check()                 # full R CMD check (use --as-cran for release work)
styler::style_pkg()               # or `air format` (air.toml present)
lintr::lint_package()             # lint
```

Scaffolding: `usethis::use_r()`, `use_test()`, `use_package()`. testthat
3e; `vdiffr` for plot snapshots; roxygen2 for every exported function
(document the *why* of each default, runnable `@examples` on `sim_dyad`,
`@seealso` cross-links).

## Definition of done (every change)

- Tests written/updated and passing; new behavior has a test.
- **C++ changed?** Recompiled via `load_all()`;
  [`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
  re-run and the regenerated `RcppExports.cpp`/`RcppExports.R` committed
  *clean* (no stray reordering churn); a **numerical-regression test**
  guards any optimization (Invariant 5); **no build artifacts** (`*.o`,
  `*.so`, `*.dll`) staged or committed.
- **Efficiency change?** `bench/` script records before/after timings;
  the milestone cites the measured speedup.
- `devtools::document()` run if roxygen changed; NAMESPACE committed.
- `devtools::check()` clean (`--as-cran` for release-track work); notes
  triaged.
- Styled (`air`/styler) and linted.
- User-visible change reflected in NEWS.md (once it exists) and the
  relevant `@examples`/vignette.

## Git

- Default branch is **`main`**.
- Small, focused commits with imperative messages (e.g.,
  `Honor na.rm in calc_wcc_cpp`).
- Don’t force-push `main`. **Don’t commit data, credentials, or build
  artifacts** (`src/*.o`, `src/*.so`, `src/*.dll`, `.DS_Store`,
  `Rplots.pdf`) — M3 untracks the ones currently slipped in.
- Commit or push only when asked.

## Ask-first / guardrails

- Ambiguity in `DESIGN.md` → ask; don’t invent a design decision.
- Adding an `Imports` dependency, changing a resolved default, or
  changing the numeric output of a C++ core (beyond the deliberate M1
  window-semantics fix) → flag for approval first.
- Adopting OpenMP / introducing C-level parallelism (M2) → flag; default
  must stay serial and reproducible.
- Touching git history, tags, or anything destructive → confirm first.
- Prefer extending the existing C++ cores and R helpers over
  reimplementing numerics or adding deps.

## Out of scope for now

- **New estimators** (phase synchrony, wavelet coherence, CRQA/MEA) —
  deferred to M6–M8, and only after the M5 shared-surface framework
  lands.
- **Group-level / multivariate modeling** — deferred to M9 (`mvSUSY` is
  the multivariate reference).
- **tidy/glance/as_tibble methods** — resolved to add (DESIGN.md
  §7/§14), built in **M5**; not before.
- **IAAFT / segment-shuffling surrogates** — resolved to add (DESIGN.md
  §6/§14), scheduled **M10**.
- **A unified `bsync_ts` preprocessing object** — logged in DESIGN.md
  §14; specify before building.
- **OpenMP** — not adopted until the M2 decision; no `#pragma omp` ships
  before then.

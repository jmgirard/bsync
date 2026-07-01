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

## Milestone bookkeeping (single source of truth)

Detailed milestone history — the baseline and the full `M<N> (done)`
narrative for every completed milestone — lives in **`MILESTONES.md`**
(repo root), in numeric order with no gaps. This file carries only the
one-line index below. To avoid redundancy, each kind of note has exactly
one home:

- **Detailed milestone narrative** → `MILESTONES.md`.
- **One-line index per milestone** → “## Completed milestones” below.
- **User-facing changes** → `NEWS.md`.
- **Design-contract changes** → `DESIGN.md` §1–14 (its §15 is the
  forward roadmap plus a pointer to `MILESTONES.md`, **not** a second
  copy of the log).

## Completed milestones

One line each; full detail in
[`MILESTONES.md`](https://jmgirard.github.io/bsync/MILESTONES.md).

- **M1 — Correctness & robustness.** `na.rm` honored in WCC; window-size
  semantics fixed (`w_max = window_size - 1`); short-series aborts;
  [`evaluate_signal_power()`](https://jmgirard.github.io/bsync/reference/evaluate_signal_power.md)
  side-effect plotting removed; cli condition style;
  [`leadership_asymmetry()`](https://jmgirard.github.io/bsync/reference/leadership_asymmetry.md)/[`suggest_wcc_params()`](https://jmgirard.github.io/bsync/reference/suggest_wcc_params.md)
  docs.
- **M2 — Efficiency.** `calc_wcc_cpp` rewritten to an NA-aware
  prefix-sum core (5.4×–24.9×); pure-R oracle regression test; OpenMP
  removed (serial, reproducible).
- **M3 — CRAN readiness.** Build artifacts untracked; `@return`
  everywhere; `R CMD check --as-cran` 0/0/0; spell-check enforced in CI.
- **M4 — Selectable WCC aggregate statistic.**
  `wcc(statistic = c("mean_abs_z","peak"))`; single `wcc_aggregate()`
  drives observed + surrogate paths (Invariant 2).
- **M5 — Shared windowed-surface + surrogate framework + tidy
  interface.** `build_surface_grid()` / `run_surrogate_engine()` /
  `build_surface_heatmap()` factored; `bsync_surface` superclass;
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html)/[`glance()`](https://generics.r-lib.org/reference/glance.html)/[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
  (`generics` + `tibble` added to Imports).
- **M6 — Parameter guidance & synchrony multiverse.**
  [`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md)
  engine + spec-curve plot;
  [`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md)
  rebuilt as a thin wrapper + exported
  [`select_specification()`](https://jmgirard.github.io/bsync/reference/select_specification.md);
  PSD-driven
  [`suggest_wcc_params()`](https://jmgirard.github.io/bsync/reference/suggest_wcc_params.md).
- **M7 — First CRAN release (`v0.1.0`).** Docs/messaging pass
  (get-started vignette, grouped pkgdown, README); version bump to
  `0.1.0`; `cran-comments.md`; ready to submit (human-gated).

## Current focus

**Post-release maintenance and next-estimator work.** M7 (`v0.1.0`) is
complete and ready to submit to CRAN (human-gated — see
`cran-comments.md` pre-submission checklist). The next milestone is M8.

- **M8 — Phase synchrony estimator (next).** Add a windowed
  phase-synchrony estimator to the `bsync_surface` framework established
  in M5. See `DESIGN.md` §15 for scope.

See `DESIGN.md` §15 for the full roadmap (M8 phase synchrony; M9 wavelet
coherence; M10 CRQA/MEA; M11 group-level workflow; M12 expanded
surrogates).

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
`future.apply`, `generics`, `ggplot2`, `grDevices`, `gsignal`, `rlang`,
`scales`, `tibble`, `utils`. `Suggests`: `future`, `knitr`, `rmarkdown`,
`spelling`, `testthat (>= 3.0.0)`, `vdiffr`. `generics` and `tibble`
added in M5 for the tidy interface. **Do not grow `Imports` without
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
- **Milestone work happens on a feature branch, merged via PR — not
  committed directly to `main`.** `/plan-milestone` cuts `m<N>-<slug>`
  off an up-to-date `main`; the planning commit (the “## Current focus”
  update) and all implementation commits land there. `main` is updated
  only by merging a reviewed PR once CI is green (see below). Trivial,
  isolated doc-typo fixes may still go directly to `main` at the user’s
  discretion; anything touching `R/`, `src/`, `tests/`, `DESCRIPTION`,
  or vignettes goes through a PR.
- **Commit / push cadence (respects “only when asked”):** the milestone
  skills create the branch and make commits **locally**; **pushing the
  branch and opening the PR happen only on the user’s request** — never
  auto-pushed. Don’t force-push `main`.
- **CI gate:** a milestone PR merges into `main` only after the three
  GitHub Actions workflows go green — `R-CMD-check`, `test-coverage`,
  `pkgdown`.
- **R-hub (on demand, not a per-PR gate):** for milestones that touch
  `src/` and before any actual CRAN submission, also run an on-demand
  R-hub check (`rhub::rhub_check()`, workflow
  `.github/workflows/rhub.yaml`) for sanitizers (ASAN/UBSAN), valgrind,
  and the extra platforms the local macOS + `R-CMD-check` matrix can’t
  see — the highest-value check for the C++/Armadillo cores. It is slow
  and `workflow_dispatch`-triggered, so it is deliberately **not** part
  of the per-PR gate.
- Small, focused commits with imperative messages (e.g.,
  `Honor na.rm in calc_wcc_cpp`).
- **Don’t commit data, credentials, or build artifacts** (`src/*.o`,
  `src/*.so`, `src/*.dll`, `.DS_Store`, `Rplots.pdf`) — M3 untracks the
  ones currently slipped in.

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
  deferred to M8–M10, and only after the M5 shared-surface framework
  lands.
- **Parameter-guidance overhaul**
  ([`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md)
  engine,
  [`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md)
  as a thin wrapper over it, PSD-driven
  [`suggest_wcc_params()`](https://jmgirard.github.io/bsync/reference/suggest_wcc_params.md))
  — deferred to **M6**, after the M5 framework (DESIGN.md §14 \#10).
- **First CRAN release** (`v0.1.0`) — explicitly its own milestone
  **M7**, after M6; not near-term.
- **Group-level / multivariate modeling** — deferred to M11 (`mvSUSY` is
  the multivariate reference).
- **tidy/glance/as_tibble methods** — built in **M5** (done).
- **IAAFT / segment-shuffling / pseudo-dyad surrogates** — resolved to
  add (DESIGN.md §6/§14), scheduled **M12**; the pseudo-dyad
  (between-dyad rMEA) generator depends on M11’s `dyad_list`.
- **A unified `bsync_ts` preprocessing object** — logged in DESIGN.md
  §14; specify before building.
- **OpenMP** — not adopted until the M2 decision; no `#pragma omp` ships
  before then.

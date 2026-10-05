# CLAUDE.md — bsync

Operating manual for AI-assisted development of this package. Read
`cairn/DESIGN.md` first and treat it as the **source of truth** for all
design decisions; this file covers *how we work*, not *what we’re
building*. When this file and `DESIGN.md` disagree, `DESIGN.md` wins for
design and this file wins for process — and flag the conflict.

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
`cairn/DESIGN.md`.

The package must serve **newcomers** (safe, loud defaults + guidance)
and **experts** (every hyperparameter exposed, raw surface accessible,
full surrogate machinery) at once.

## Milestone history (pre-cairn)

Milestones M1–M7 (through the v0.1.0 CRAN-ready release) were tracked in
a precursor system; their full narrative is entombed verbatim in
`cairn/legacy/MILESTONES.md`. Status and the forward roadmap now live in
`cairn/ROADMAP.md` (see “Project tracking (cairn)” below).

**User-facing docs must not reference milestone numbers.** README,
vignettes, roxygen (and the generated `man/*.Rd`), and `NEWS.md`
describe features in plain terms — never `M<n>`. `NEWS.md` is organized
by release/theme, not by milestone. Internal `#` code comments in `R/`
*may* carry milestone tags for dev provenance.
`tests/testthat/test-doc-hygiene.R` enforces this (scans the user-facing
docs, skips `R/`).

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

## Resolved defaults (see `cairn/DESIGN.md` §9)

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

## Dependencies (see `cairn/DESIGN.md` §10)

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
- Branching, commits, merge approval, and push cadence follow the cairn
  git model (see “Project tracking (cairn)” below): milestone work on
  `m<nnn>-<slug>` branches merged via PR at the review gate, never
  directly to `main`.
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

- Ambiguity in `cairn/DESIGN.md` → ask; don’t invent a design decision.
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
  — deferred to **M6**, after the M5 framework (cairn/DESIGN.md §14
  \#10).
- **First CRAN release** (`v0.1.0`) — explicitly its own milestone
  **M7**, after M6; not near-term.
- **Group-level / multivariate modeling** — deferred to M11 (`mvSUSY` is
  the multivariate reference).
- **tidy/glance/as_tibble methods** — built in **M5** (done).
- **IAAFT / segment-shuffling surrogates** — resolved to add
  (cairn/DESIGN.md §6/§14); a candidate row in `cairn/ROADMAP.md`. The
  pseudo-dyad (between-dyad rMEA) generators shipped in M009 on the
  existing `dyad_list` input.
- **A unified `bsync_ts` preprocessing object** — logged in
  cairn/DESIGN.md §14; specify before building.
- **OpenMP** — not adopted until the M2 decision; no `#pragma omp` ships
  before then.

## Project tracking (cairn)

This repo uses the cairn plugin. **Before acting on any request,
classify it and route** — the tracking rulebook only loads once a cairn
skill fires, so starting work in plain conversation silently bypasses
the work tiers and the git model. Classify first:

- **Trivial** (no runtime surface — typo, comment, tracking edit):
  commit directly to the default branch.
- **User-visible bug**: invoke `/hotfix`.
- **New work, a design decision, or more than one sitting**: invoke
  `/milestone-plan` (then `/milestone-implement` → `/milestone-review`).
- **Status, “what’s next”, or unsure which tier**: invoke `/milestone`.
- **Never implement code on the default branch** outside a
  milestone/hotfix branch; nothing reaches it without the user’s
  explicit approval at the review gate.

Whenever the request is anything but trivial, invoke the skill *first*
so the full rulebook (the plugin’s `skills/shared/tracking-rules.md`)
and its conduct load — do not reconstruct the rules here from memory.
All project state lives under `cairn/` (**Architecture → DESIGN · Status
→ ROADMAP · Tasks → milestone files · Decisions → DECISIONS · Lessons →
LESSONS · History → archive + git**); never record status or TODOs in
this file. Claude’s persistent memory never holds project state;
`cairn/` files win any conflict.

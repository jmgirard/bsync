# CLAUDE.md — `bsync`

Operating manual for AI-assisted development of this package. Read `DESIGN.md` (repo root) first and
treat it as the **source of truth** for all design decisions; this file covers *how we work*, not
*what we're building*. When this file and `DESIGN.md` disagree, `DESIGN.md` wins for design and this
file wins for process — and flag the conflict.

## What this is

`bsync` is an R package for analyzing **interpersonal / behavioral synchrony** in continuous dyadic
time series. It provides C++/Armadillo-accelerated windowed estimators of nonstationary lead–lag
structure (windowed cross-correlation, windowed dynamic time warping, windowed Granger causality), a
preprocessing pipeline (PSD-based downsampling guidance, zero-phase smoothing, kinematics, gap
imputation, time-bin aggregation), surrogate (pseudo-synchrony) significance testing, peak/valley
optima picking and leadership-asymmetry indices, and theory- and data-driven hyperparameter helpers.
Full rationale, contracts, the surface object spec, and resolved defaults are in `DESIGN.md`.

The package must serve **newcomers** (safe, loud defaults + guidance) and **experts** (every
hyperparameter exposed, raw surface accessible, full surrogate machinery) at once.

## Baseline (pre-milestone)

Substantial functionality predates formal milestone tracking; new milestones begin at **M1**. The
baseline, as built:

- **Estimators:** `wcc()`, `wdtw()`, `wgranger()` with Rcpp/Armadillo cores
  (`calc_wcc_cpp`, `calc_wdtw_cpp`, `calc_wgranger_cpp`) and S3 `print`/`summary`/`plot`.
- **Optima & leadership:** `pick_optima()` (+ `pick_optima_cpp`), `leadership_asymmetry()`,
  `plot_optima_overlay()`.
- **Surrogates:** `generate_surrogate_circular()`, `generate_surrogate_phase()`,
  `wcc_surrogate()`, `wdtw_surrogate()`, `wgranger_surrogate()` (parallel via `future.apply`).
- **Preprocessing:** `evaluate_signal_power()`, `downsample_signal()`, `aggregate_by_time()`,
  `smooth_signal()`, `trim_edges()`, `calc_velocity_1d()`, `calc_speed_{1,2,3}d()`,
  `diagnose_ts_gaps()`, `impute_ts_gaps()`.
- **Tuning:** `suggest_wcc_params()`, `autotune_wcc()`.
- **Data:** `sim_dyad` (simulated 3D dyadic positions with a shifting lead–lag relationship).
- **Infra:** testthat 3e suite, `vdiffr` plot snapshots, pkgdown config, vignettes (WCC workflow,
  WDTW, Granger, surrogate testing, downsampling, WCC params).

A read of the baseline surfaced the defects M1–M3 address (see Current focus and `DESIGN.md` §14/§15).

## Completed milestones

- **M1 — Correctness & robustness (done).** All six acceptance criteria met; post-review fixes
  applied (commits `741bf7c`–`ed5960f` on `main`; 349 tests passing, 0 errors/0 warnings/2 notes
  in `R CMD check`; notes are M3 scope):
  1. `na.rm` honored in WCC: `calc_wcc_cpp` gains `na_rm` bool; `na.rm = FALSE` returns `NA` for
     any window containing `NA`; forwarded through `create_wcc_df()` and `wcc_surrogate()`. Both
     modes tested in `test-wcc.R`.
  2. Window-size semantics fixed to exactly `window_size` samples via `w_max = window_size - 1` at
     the C++ boundary in all three estimators and all three surrogate wrappers; `n_r` grid math
     updated; a test asserts realized window count (33, not 32, for the reference parameters).
  3. Short-series robustness: `seq_len()` throughout; `cli::cli_abort()` in all three
     `create_*_df()` builders and all three surrogate grid builders; abort tested for all three
     estimators.
  4. `evaluate_signal_power()` no longer auto-prints: `plot=` arg removed; `plot.signal_power_res()`
     S3 method added; no `Rplots.pdf` side effect; other compute functions audited.
  5. Condition style unified to `cli` in `R/impute.R` and `R/surrogate_generation.R`.
  6. `leadership_asymmetry()` roxygen states centered sliding-window semantics; `min_valid`
     validated; `suggest_wcc_params()` 4-cycles-per-window heuristic documented in `@details`.
  Post-review: NEWS.md header fixed; OpenMP flags removed from Makevars pending M2 decision;
  surrogate `@details` document null-statistic semantics and `fast_method` caveat; surrogate
  robustness and p-value sanity tests added (349 total).

- **M2 — Efficiency (done).** All five acceptance criteria met; post-review fixes applied
  (commits `ce385b1`–`644d140` on `main`; 353 tests passing, 0 errors/0 warnings/0 notes in
  `R CMD check`):
  1. `calc_wcc_cpp` rewritten to an NA-aware prefix-sum algorithm (loops over distinct lags,
     builds six masked prefix arrays in O(n) per τ, evaluates each window in O(1)). Both `na_rm`
     modes, bounds-check→NA, and variance-zero guard preserved. Signature unchanged; `RcppExports`
     regenerates with no diff.
  2. Pure-R `stats::cor` oracle matches new core at 1e-9 on `sim_dyad` for clean, NA na.rm=TRUE,
     and NA na.rm=FALSE cases (`test-wcc.R`). Oracle was first validated against the pre-rewrite
     core, proving it correct before the optimization.
  3. `bench/bench_wcc.R` + `bench/RESULTS.md` record 5.4×–24.9× speedup across four configs
     (sim_dyad narrow/wide, n=10000 narrow/wide). Speedup grows with w_max as expected.
  4. OpenMP removed: no `SHLIB_OPENMP` flags in `Makevars`/`Makevars.win`, no `#pragma omp`
     in any source file. Decision logged in `DESIGN.md §14`; core is serial and reproducible.
  5. 353 tests green; vdiffr snapshots unchanged; no build artifacts staged.
  Post-review: explicit C++ stdlib includes added (`<algorithm>`, `<unordered_map>`, etc.) for
  CRAN portability; large-mean oracle test added documenting ~2e-6 prefix-sum cancellation loss
  (tolerance 1e-5 catches real bugs); `.Rbuildignore` gains `^\.claude$`, `^CLAUDE\.md$`,
  `^DESIGN\.md$`, `^bench$` — `R CMD check` is now 0/0/0.

- **M3 — CRAN readiness (done).** All seven acceptance criteria met (commits `184e874`–`69f1ca4`
  on `main`; 353 tests passing, 0 errors/0 warnings/0 notes in `R CMD check --as-cran`). No C++
  core numerics changed — only `RcppExports` regenerated to confirm zero diff, so Invariant 5 was
  not triggered. Plan-time decisions held: version stays `0.0.0.9000`; `cran-comments.md` deferred
  to actual submission; DESCRIPTION Description refreshed:
  1. Build artifacts untracked via `git rm --cached` (`.DS_Store`, `tests/.DS_Store`, six
     `src/*.o`, `src/bsync.{so,dll}`, `tests/testthat/Rplots.pdf`); `.gitignore` extended with
     `*.o`/`*.so`/`*.dll`/`*.dylib` + `tests/testthat/Rplots.pdf`. `git ls-files` is artifact-free;
     `bsync.Rproj`/`LICENSE.md` deliberately kept (already `.Rbuildignore`'d).
  2. `Rcpp::compileAttributes()` produced zero diff in `src/RcppExports.cpp` / `R/RcppExports.R`.
  3. `@return` roxygen added to all 18 exported `print`/`plot`/`summary` methods +
     `plot_optima_overlay`; re-documented. Only `sim_dyad` (data) and `bsync-package` (overview)
     omit `\value`, both legitimately.
  4. `R CMD check --as-cran` = 0/0/0 (verified twice during implementation); examples needed no
     `\donttest{}` (check reports "examples ... NONE").
  5. `urlchecker::url_check()` clean; `spelling::spell_check_package()` clean via new `inst/WORDLIST`
     and `Language: en-US` in DESCRIPTION (neither tool added to deps).
  6. `devtools::build_readme()` regenerated `README.md` (and fixed a broken `wcc_surrogate()`
     example that passed a non-existent `n_surrogates` arg instead of a `y_surrogates` matrix);
     `pkgdown::check_pkgdown()` passes; DESCRIPTION Description now covers WDTW/Granger/surrogates/
     optima/leadership.
  7. 353 tests green; vdiffr snapshots unchanged; styler applied across tests + vignettes; NEWS.md
     gained an M3 entry.
  Post-review: deprecated `context()` removed from `test-impute.R` (suite now 0 warnings);
  spell-check enforced in CI via `usethis::use_spell_check()` (adds `tests/spelling.R` skip-on-cran
  + `spelling` to `Suggests`) so `inst/WORDLIST` can no longer drift unnoticed — reverses the
  plan-time "ad hoc, not a dep" choice by design. `R CMD check --as-cran` remains 0/0/0.

- **M4 — Selectable WCC aggregate statistic (done).** All seven acceptance criteria met (commits
  `0eafe21`–`7352eac` on `main`; 363 tests passing, 0 errors/0 warnings/0 notes in
  `R CMD check --as-cran`). No C++ core change — Invariants 5/6 not triggered; `RcppExports` diff
  is empty:
  1. `wcc(statistic = c("mean_abs_z","peak"))` with `match.arg`; default `"mean_abs_z"` reproduces
     pre-M4 `fisher_z` value bit-for-bit on `sim_dyad`; invalid value aborts. Tested.
  2. `peak` = mean-over-windows of max-|Fisher-z|-across-lags, verified against an independent
     pure-R oracle on `sim_dyad` to 1e-9. Tested.
  3. Single `wcc_aggregate(z, window_id, statistic)` helper (refactored from `fisher_z()`) drives
     both observed and surrogate paths; `wcc_surrogate(statistic="peak")$observed_z` equals
     `wcc(..., statistic="peak")$fisher_z` exactly — matched null (Invariant 2); tested for both
     statistics.
  4. `print.wcc_res` / `print.wcc_surr` label the aggregate by chosen statistic.
  5. `vignettes/wcc-workflow.Rmd` gains §3.1 documenting `mean_abs_z` (SUSY) vs `peak` (rMEA/
     Boker) with lineage citations and a code example. NEWS.md gained an M4 entry. `inst/WORDLIST`
     updated (rMEA, SUSY, Selectable); CI spell-check stays green.
  6. No C++ change; no new `Imports`; no build artifacts staged.
  7. 363 tests green; `R CMD check --as-cran` = 0/0/0; styled (styler) and linted.
  Plan-time decisions held: WCC-only scope; `$fisher_z` field name kept (revisit in M5);
  `settings$statistic` records the choice.
  Post-review: peak oracle rewritten to an explicit for-loop (no `tapply`) so it cannot
  share a latent bug with the implementation; cross-path Invariant-2 test added (pass `y`
  itself as the sole surrogate, assert `surrogate_z[1] == observed_z` for both statistics —
  exercises the surrogate loop's `wcc_aggregate(grid_df$row)` call site independently);
  peak + all-NA-window, peak + `time`-supplied grouping, peak print-label, and peak
  p-value sanity tests added (373 tests total, `R CMD check --as-cran` remains 0/0/0).

## Current focus

**Hardening + core-completion cycle toward a first CRAN release**, run as focused milestones via the
plan → implement → review loop. The first release is **explicitly not near-term** — it is its own
milestone (M7, `v0.1.0`) after the M5 framework and M6 parameter guidance land, not a deadline
hanging over the hardening work. **M5 is the active milestone.**

- **M5 — Shared windowed-surface + surrogate framework + tidy interface (active).** Factor the
  duplicated grid-builder, validation, surrogate engine, and plot scaffold into one shared layer so
  new estimators plug in cheaply, and add the broom-style tidy interface. **No C++ change**
  (RcppExports diff must be empty); the refactor is output-preserving (the no-C++ analogue of
  Invariant 5) — observed numerics, `results_df` schema, `print`/`summary` text, and all vdiffr
  snapshots stay byte-identical.

  *Resolved at plan time:* (1) **Granger → shared machinery, estimator-specific surface.** All three
  estimators share the validator, grid builder, surrogate engine, `bsync_surface` superclass, and
  tidy layer, but keep their own `results_df` columns (WCC/WDTW lagged with `tau`; Granger lag-free,
  directional `f_xy/p_xy/f_yx/p_yx`) and their own plot geometry (WCC/WDTW share the heatmap
  scaffold; Granger keeps its line plot). Granger is *not* contorted into a symmetric similarity
  (resolves DESIGN.md §14 "remaining"). (2) **`bsync_surface` superclass** — `wcc_res`/`wdtw_res`/
  `wgranger_res` inherit it; shared/tidy methods dispatch on it, estimator-specifics on the leaf.
  (3) **Uniform `$aggregate` slot** — drop `$fisher_z`/`$mean_distance` (no users, no shim); all
  surface objects carry a named-numeric `$aggregate` (`mean_abs_z`/`mean_distance`/`c(f_xy, f_yx)`).
  (4) **`generics` + `tibble` → Imports** (tibble already transitive via dplyr; no new footprint).

  *Acceptance criteria:*
  1. **One grid builder.** `build_surface_grid()` (in `R/surface.R`) is the sole site of the `n_r`
     math, the `w_max = window_size - 1` boundary (Inv 4), the lag sequence, and the short-series
     abort; all three `create_*_df` builders **and** all three surrogate wrappers call it (the
     surrogate wrappers currently re-derive the grid a fourth time — that copy dies). `grep` finds
     the `floor((n_x - w_max ...))` arithmetic in exactly one place. `wcc()`/`wdtw()`/`wgranger()`
     `results_df` + aggregates are bit-identical to pre-refactor on `sim_dyad` (characterization
     test captured before the refactor).
  2. **One validator.** A shared validation path replaces the per-estimator x/y/time/integerish
     blocks; abort conditions/messages preserved and tested.
  3. **`bsync_surface` superclass + uniform `$aggregate`.** All three objects inherit
     `"bsync_surface"` and carry `$aggregate`; existing `print`/`summary` text and all vdiffr
     snapshots unchanged.
  4. **One surrogate engine.** `run_surrogate_engine()` drives all three wrappers; observed and
     surrogate statistics flow through the same grid builder + aggregate fn (Inv 2/4); tail
     direction parameterized once (upper WCC/Granger, lower WDTW). The M4 cross-path Inv-2 test
     (pass `y` as its sole surrogate ⇒ null draw == observed) is extended to WDTW and Granger (both
     directions). p-values bit-identical to pre-refactor for a fixed seed. **M6-enabling seams
     (contractual, not incidental):** the engine accepts a **prebuilt grid** and a **prebuilt
     surrogate matrix**, and exposes an **aggregate-only path** (core → aggregate, no per-cell
     `results_df`), so the M6 multiverse can hoist surrogate generation and grid construction out
     of its inner loop and reuse one surrogate matrix across every parameter cell sharing a
     `surrogate_method`. A test asserts the engine never materializes a full `results_df` on the
     surrogate path.
  5. **One heatmap scaffold.** `build_surface_heatmap()` factors the `time_step`/axis/theme/
     zero-line scaffold; `plot.wcc_res`/`plot.wdtw_res` supply only the fill scale; all plot
     snapshots unchanged. Granger retains its line plot.
  6. **Tidy interface.** `tidy()`/`glance()`/`as_tibble()` registered on `bsync_surface` for all
     three estimators; `tidy()` = long surface, `glance()` = one-row run meta (estimator,
     window/lag/ar settings, aggregate(s), n_windows[, n_lags]), `as_tibble()` ≡ `tidy()`; Granger's
     tau-less / two-aggregate shape handled. `generics` + `tibble` in Imports; NAMESPACE
     regenerated. Tested for all three.
  7. **Docs & DESIGN.** NEWS.md M5 entry; tidy-interface vignette/section on `sim_dyad`; WORDLIST
     updated; DESIGN.md §14 Granger/superclass item moved from "remaining" to "resolved".
  8. **Gates.** Full suite green (current 373 + new); `R CMD check --as-cran` 0/0/0; styler + lintr
     clean; **no C++ change** (RcppExports regenerated, diff empty); no build artifacts staged.

- **M6 — Parameter guidance & synchrony multiverse.** One engine, three read-outs on the M5 shared
  framework: `synchrony_multiverse()` is the grid + matched-null-surrogate engine (headline metric =
  ES vs. null, not raw synchrony) with a specification-curve plot; `autotune_wcc()` becomes a thin
  wrapper = multiverse + a selection rule (detectability + cross-dyad stability, not bare ES argmax),
  validated against `sim_dyad`; `suggest_wcc_params()` stays the single PSD-data-driven starting
  point with the SUSY constraints enforced/reported. There is no single "correct" WCC parameter set —
  the optimum depends on the signal's own timescales — so the honest deliverable is the multiverse,
  with the matched-null surrogate as the defense against autocorrelation-driven spurious correlation
  (DESIGN.md §14 #10, §15 M6). Engine grid axes are in **time units (seconds), not samples**, so
  downsample/smoothing can be *opt-in* multiverse axes (default fixes preprocessing). The
  `autotune_wcc()` validation lands here, before the M7 release.

- **M7 — First CRAN release (`v0.1.0`).** Cut the first public release once WCC/WDTW/Granger + the M5
  framework + tidy interface + M6 parameter guidance cohere; later estimators are post-1.0 minors.

See `DESIGN.md` §15 for the full roadmap (M7 first CRAN release; M8 phase synchrony; M9 wavelet
coherence; M10 CRQA/MEA; M11 group-level workflow; M12 expanded surrogates).

## Invariants — do not violate without flagging

These encode hard-won reasoning. Changing them is a design decision, not a refactor.

1. **C++ cores are pure and validated in R.** All argument checking, NA policy, and grid
   construction live in the R wrapper; cores assume clean inputs, do their own bounds checks
   (returning `NA` out of range), and never message the user. Never expose a `*_cpp` function as
   user API.
2. **Surrogate nulls match the observed statistic.** The aggregate statistic computed on the
   observed data and on every surrogate must be identical; the p-value is its tail. Change one,
   change both in lockstep.
3. **Loud, non-destructive preprocessing.** Diagnostics advise; action functions transform only what
   is asked and never silently change the basis. Announce consequential auto-choices via cli.
4. **`window_size` means exactly `window_size` samples.** Realized window length is a contract;
   `w_max = window_size - 1` at the C++ boundary is its single source of truth.
5. **Optimizations never change results.** Any change to a C++ core reproduces the prior
   implementation's output within tolerance on `sim_dyad` (and an NA-laden case). Speed is free;
   numbers are sacred.
6. **Reproducible stochastics.** Surrogate generation and autotune sampling respect `set.seed()` and
   `future.seed = TRUE`; never reseed internally.
7. **Light result objects.** Carry `results_df` + `settings` + the aggregate; raw surrogate draws
   live only in surrogate objects; raw input data is not stored.
8. **Time integrity.** When `time` is supplied, window positions map to real timestamps; edge
   trimming preserves the true timeline (the documented reason `time` exists).

## Resolved defaults (see `DESIGN.md` §9)

`window_size`/`lag_max` **required** · window length = exactly `window_size` samples · increments =
`1`/`1` · WCC `na.rm = TRUE` (pairwise; honored — `FALSE` ⇒ NA window) · WDTW `scale_method =
"global"`, `distance_metric = "L2"` · Granger `ar_order = 1` · WCC aggregate statistic selectable,
default `"mean_abs_z"` (M4) · surrogate method user-chosen (phase = spectrum, circular =
autocorrelation) · `n_surrogates = 100` (≥ 1000 advised for reporting) · smoothing = Savitzky–Golay
order 3 · downsample/aggregate = median · `impute maxgap = 5`, no extrapolation · PSD `threshold =
0.95`. Don't change these silently.

## Dependencies (see `DESIGN.md` §10)

Compiled cores via **`LinkingTo: Rcpp, RcppArmadillo`** are the package's reason for existing —
unlike a pure-R package, heavy inner loops belong in C++ here. Current `Imports`: `Rcpp`, `cli`,
`dplyr`, `future.apply`, `ggplot2`, `grDevices`, `gsignal`, `rlang`, `scales`, `utils`. `Suggests`:
`future`, `knitr`, `rmarkdown`, `testthat (>= 3.0.0)`, `vdiffr`. **Do not grow `Imports` without
flagging it** — prefer base R or the existing stack, and put heavy compute in the C++ cores, not new
R deps.

## Dev workflow

R ≥ 4.1 (native pipe `|>` and `\(x)` lambdas allowed). Standard devtools loop, **with the compiled
step**:

```r
devtools::load_all()              # compiles changed C++ + loads; run after any src/ edit
Rcpp::compileAttributes()         # regenerate RcppExports after changing a [[Rcpp::export]] signature
devtools::document()              # regenerate roxygen docs + NAMESPACE after any roxygen change
devtools::test()                  # run testthat suite
devtools::check()                 # full R CMD check (use --as-cran for release work)
styler::style_pkg()               # or `air format` (air.toml present)
lintr::lint_package()             # lint
```

Scaffolding: `usethis::use_r()`, `use_test()`, `use_package()`. testthat 3e; `vdiffr` for plot
snapshots; roxygen2 for every exported function (document the *why* of each default, runnable
`@examples` on `sim_dyad`, `@seealso` cross-links).

## Definition of done (every change)

- Tests written/updated and passing; new behavior has a test.
- **C++ changed?** Recompiled via `load_all()`; `Rcpp::compileAttributes()` re-run and the
  regenerated `RcppExports.cpp`/`RcppExports.R` committed *clean* (no stray reordering churn); a
  **numerical-regression test** guards any optimization (Invariant 5); **no build artifacts**
  (`*.o`, `*.so`, `*.dll`) staged or committed.
- **Efficiency change?** `bench/` script records before/after timings; the milestone cites the
  measured speedup.
- `devtools::document()` run if roxygen changed; NAMESPACE committed.
- `devtools::check()` clean (`--as-cran` for release-track work); notes triaged.
- Styled (`air`/styler) and linted.
- User-visible change reflected in NEWS.md (once it exists) and the relevant `@examples`/vignette.

## Git

- Default branch is **`main`**.
- Small, focused commits with imperative messages (e.g., `Honor na.rm in calc_wcc_cpp`).
- Don't force-push `main`. **Don't commit data, credentials, or build artifacts** (`src/*.o`,
  `src/*.so`, `src/*.dll`, `.DS_Store`, `Rplots.pdf`) — M3 untracks the ones currently slipped in.
- Commit or push only when asked.

## Ask-first / guardrails

- Ambiguity in `DESIGN.md` → ask; don't invent a design decision.
- Adding an `Imports` dependency, changing a resolved default, or changing the numeric output of a
  C++ core (beyond the deliberate M1 window-semantics fix) → flag for approval first.
- Adopting OpenMP / introducing C-level parallelism (M2) → flag; default must stay serial and
  reproducible.
- Touching git history, tags, or anything destructive → confirm first.
- Prefer extending the existing C++ cores and R helpers over reimplementing numerics or adding deps.

## Out of scope for now

- **New estimators** (phase synchrony, wavelet coherence, CRQA/MEA) — deferred to M8–M10, and only
  after the M5 shared-surface framework lands.
- **Parameter-guidance overhaul** (`synchrony_multiverse()` engine, `autotune_wcc()` as a thin
  wrapper over it, PSD-driven `suggest_wcc_params()`) — deferred to **M6**, after the M5 framework
  (DESIGN.md §14 #10).
- **First CRAN release** (`v0.1.0`) — explicitly its own milestone **M7**, after M6; not near-term.
- **Group-level / multivariate modeling** — deferred to M11 (`mvSUSY` is the multivariate reference).
- **tidy/glance/as_tibble methods** — resolved to add (DESIGN.md §7/§14), built in **M5**; not before.
- **IAAFT / segment-shuffling / pseudo-dyad surrogates** — resolved to add (DESIGN.md §6/§14),
  scheduled **M12**; the pseudo-dyad (between-dyad rMEA) generator depends on M11's `dyad_list`.
- **A unified `bsync_ts` preprocessing object** — logged in DESIGN.md §14; specify before building.
- **OpenMP** — not adopted until the M2 decision; no `#pragma omp` ships before then.

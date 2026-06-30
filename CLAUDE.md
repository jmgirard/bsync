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

- **M5 — Shared windowed-surface + surrogate framework + tidy interface (done).** All nine
  acceptance criteria met (commits `0769274`–`b17b464` on `main`; 449 tests passing, 0 errors/0
  warnings/0 notes in `R CMD check --as-cran`). No C++ change — Invariants 5/6 not triggered;
  `RcppExports` diff is empty:
  0. External-oracle preflight: WDTW L1 golden (6.22…, dtw v1.23-3 symmetric1) and Granger F
     golden (f_xy=287.62…, lmtest v0.9-40) committed as frozen fixtures in `test-external-oracle.R`;
     WCC oracle uses pre-existing base-`cor` path. All 8 preflight tests pass.
  1. `build_surface_grid()` is the sole site of `n_r` math and `w_max = window_size − 1`; all three
     `create_*_df` builders and all three surrogate wrappers call it. `wcc()`/`wdtw()`/`wgranger()`
     `results_df` + aggregates bit-identical to pre-refactor on `sim_dyad` (characterization freeze
     after preflight). Validated via 41 AC1/AC2/AC3 characterization and validator tests.
  2. `validate_series()` + `validate_window_params()` replace per-estimator x/y/time/integerish
     blocks; abort conditions preserved and tested.
  3. All three objects inherit `"bsync_surface"` and carry a named-numeric `$aggregate`;
     `$fisher_z`/`$mean_distance` dropped (no users, no shim); all `print`/`summary` text and
     vdiffr snapshots unchanged.
  4. `run_surrogate_engine()` drives all three wrappers via aggregate-only closures (no `results_df`
     on surrogate path — Invariant 7 enforced); accepts prebuilt grid + surrogate matrix (M6 seam);
     cross-path Invariant-2 test extended to WDTW and Granger (both directions).
  5. `build_surface_heatmap()` factors `time_step`/axis/theme/zero-line logic;
     `plot.wcc_res`/`plot.wdtw_res` supply only the fill scale; all plot snapshots unchanged.
  6. `tidy()`/`glance()`/`as_tibble()` on `bsync_surface`: `tidy()` = long surface tibble;
     `glance()` = 1-row summary (aggregate(s) + settings + n_windows); Granger's two-aggregate
     shape handled as separate columns. `generics` + `tibble` added to `Imports`.
  7. NEWS.md M5 entry; wcc-workflow.Rmd §8 tidy-interface demo; WORDLIST updated; vignette
     `$fisher_z` references corrected. DESIGN.md §14 Granger/superclass item marked resolved.
  8. 449 tests pass (0 skip on CRAN = vdiffr); `R CMD check --as-cran` 0/0/0; styled (styler);
     no C++ change; no build artifacts staged.
  Plan-time decisions held: Granger kept estimator-specific surface (no contortion into symmetric
  similarity); `$fisher_z`/`$mean_distance` dropped rather than shimmed; `generics` + `tibble`
  added to Imports as approved.
  Post-review (commit `ecc6cf9`): four findings fixed. (1) `wdtw_surrogate(fast_method = TRUE)`
  window-grid regression — an interim refactor used a lag-free grid + `lag_max` shift that
  over-counted windows past the series end (spurious out-of-range `NA`s) and changed the fast-path
  p-value; restored to the lagged grid's distinct window starts at `tau = 0` so surrogates cover
  exactly the observed surface's windows (slow/default path was never affected; NEWS.md gains a
  bug-fix entry). (2) `glance()` used bare `%||%` (base R >= 4.4 only) while DESCRIPTION declares
  `R (>= 4.1)` — added `@importFrom rlang %||%`. (3) DESIGN.md §14 shared-surface/superclass item
  moved from "Remaining" to "Resolved (M5)" (the AC7 doc move that was claimed but not done). (4)
  AC4 test gaps closed: Invariant-7 no-`results_df` assertion on all three surrogate objects, a
  seeded p-value regression guard (WCC 91/99, WDTW 83/99 — non-boundary, exercises tail counting),
  and a `fast_method` window-alignment regression test. 455 tests pass; `R CMD check --as-cran`
  remains 0/0/0; no C++ change.

- **M7 — First CRAN release (`v0.1.0`) (done).** All nine acceptance criteria met (commits
  `633d58d`–`dcd9ef0` on `main`, plus post-review polish through `7d848c6`; 564 tests passing
  pre-polish, 0 errors/0 warnings/0 notes in `R CMD check --as-cran`). No C++ change — Invariants
  5/6 not triggered; `RcppExports` diff empty. Run as two phases under M7 number:
  **Phase A** (docs/messaging): 6 ACs met.
  1. `_pkgdown.yml` restructured with grouped `articles:` (3 groups) and `reference:` (8 groups);
     `pkgdown::check_pkgdown()` clean.
  2. New `vignettes/bsync.Rmd` "Get started" — end-to-end `sim_dyad` arc, WCC/WDTW/WGC decision
     table, reading map into 6 deep-dives.
  3. `choosing-parameters.Rmd` retitled "Choosing Analysis Parameters"; all multiverse/autotune chunks
     now run on `sim_dyad` (`n_surrogates = 100L`); `select_specification()` taught; `wgranger`
     estimator example added; `tidy()`/`glance()` shown.
  4. README regenerated: Granger added to headline; "Where to go next" article table; param-guidance
     pointer; output regenerated. Bug fix: `print.wcc_res`/`print.wcc_surr` no longer display a
     spurious double colon in the aggregate label; label renamed `Mean |Fisher's Z|` → `Mean Abs.
     Fisher's Z` to avoid cli glue-in-key-name interpolation.
  5. See-also footers added to `wdtw-workflow`, `wgranger-workflow`, `determine-downsampling`.
  6. No new `Imports`; `inst/WORDLIST` updated; spell-check green; `R CMD check --as-cran` 0/0/0.
  Post-plan deviation: AC6 stated "no R/ change" but the cli double-colon bugfix (a pre-existing
  display defect since M4) required touching `R/wcc.R` and `R/surrogate_analysis.R` — purely
  cosmetic print-method fix, no numeric/behavioral change, no C++ touched, no NAMESPACE change.
  **Phase B** (release cut): 3 ACs met.
  7. Version bump `0.0.0.9000` → `0.1.0`; NEWS.md header → `# bsync 0.1.0`; lifecycle badge
     promoted `experimental` → `stable` in README.
  8. `cran-comments.md` written (0/0/0 local; cross-platform win-builder/R-hub pre-submission
     checklist). `spelling` + `pkgdown::check_pkgdown()` clean. cran-comments.md added to
     `.Rbuildignore`. (The pkgdown site has since been deployed: `articles/bsync.html` is live and
     `urlchecker::url_check()` reports all URLs correct — the original "two 404s pending deploy"
     caveat is resolved; `cran-comments.md` updated to reflect it.)
  9. Plan stops at "ready to submit"; `print.wcc_surr` label regression test added.
  Actual CRAN upload is the user's action; pre-submission checklist in `cran-comments.md`.
  Post-review polish (commits `5c4a97c`–`7d848c6`, 564 tests, `R CMD check --as-cran` remains
  0/0/0): (1) `autotune_wcc()` now returns a classed `bsync_autotune` object with a tidy `print()`
  method (was an unclassed list that dumped the per-dyad multiverse list to the console);
  regression test added. (2) pkgdown navbar `intro` component restored so the Get-started
  `vignette("bsync")` link reappears; the "Get started" articles group given a navbar heading. (3)
  `determine-downsampling.Rmd` corrected to call `plot()` on the `evaluate_signal_power()` result
  (the stale `$plot` list access rendered no figure); `wdtw_surrogate()` `@details` reference to
  the M5-removed `$mean_distance` field corrected to `$aggregate[["mean_distance"]]`. (4) vignettes
  call the re-exported `tidy()`/`glance()` verbs directly instead of the `generics::` prefix. (5)
  root-level `Rplots.pdf` added to `.gitignore`.

- **M6 — Parameter guidance & synchrony multiverse (done).** All nine acceptance criteria met
  (commits `026b2f4`–`bb5f6bc` on `main`; 542 tests passing, 0 errors/0 warnings/0 notes in
  `R CMD check --as-cran`). No C++ change — Invariants 5/6 not triggered; `RcppExports` diff empty.
  Run as two phases under M6 number; no new `Imports`:
  1. `synchrony_multiverse(x, y, estimator, sample_rate, window_sec, lag_sec, ...)` sweeps a
     seconds-specified parameter grid (converted to samples per cell; `lag_max` hard-capped at
     `floor(window_size/2)`), runs matched-null surrogate test per cell, and returns a light
     `bsync_multiverse` (tidy grid + settings + robustness summary; Invariant 7). All three
     estimators supported (WCC/WDTW/Granger); Granger has two ES/p columns (es_xy/es_yx). Tested.
  2. ES polarity correct: WCC/Granger upper-tail `(obs-null_mean)/null_sd`; WDTW lower-tail
     `(null_mean-obs)/null_sd`. Cross-path Invariant-2 tests for all three estimators.
  3. Surrogate-reuse seam: one `y_surrogates` matrix per `surrogate_method` reused across all
     cells — asserted by spy tests (`local_mocked_bindings`); per-method cost never multiplies
     by grid size.
  4. Granger no-lag-axis and WDTW smoke + correctness tested.
  5. `print`/`summary`/`tidy`/`glance`/`as_tibble`/`plot` on `bsync_multiverse`; `plot` is a
     Simonsohn-style spec curve (ES panel + choice dashboard; pure ggplot2 + base `grid`); vdiffr
     snapshot added.
  6. `suggest_wcc_params(x, y, sample_rate, ...)` reworked: derives dominant timescale from PSD
     via `evaluate_signal_power()`; `event_duration_sec` optional override; three hard constraints
     enforced and warned (SUSY lag cap, series-length ceiling, min-samples floor). Signature change
     documented in NEWS.
  7. `autotune_wcc(dyad_list, ...)` rewritten as thin wrapper over `synchrony_multiverse()` +
     `select_specification()` (gated stability-penalized score: detectability gate + `median(ES) -
     iqr_penalty*IQR(ES)`); `select_specification()` exported for advanced users.
  8. `sim_dyad` regression test (`autotune_wcc()` on 3-dyad list, median_es > 0) and synthetic
     heterogeneous-dyad stability test both pass. Gate-fallback warning tested explicitly.
  9. 542 tests pass; `R CMD check --as-cran` 0/0/0; styled (styler); no new `Imports`; no C++
     change; all heavy tests `skip_on_cran`; `vignettes/choosing-parameters.Rmd` builds;
     NEWS.md/WORDLIST/DESIGN.md §14 #10 updated.
  Post-plan deviations: (a) `stability_flag`/`top-k` from AC7 simplified to `sig_rate`, `iqr_es`,
  `n_cells_gated` — same information, cleaner API; (b) uncoupled negative controls from AC8
  replaced by heterogeneous-noise synthetic dyads (positive-signal controls with varying SNR).
  Non-ASCII characters in R code files (em-dash, en-dash, math symbols from roxygen) were the main
  R CMD check hurdle; `%` in multiline Rd `\item` content also required `\%` escaping.
  Post-review fixes (550 tests; `R CMD check --as-cran` remains 0/0/0): (1) AC2 cross-path
  Invariant-2 gap closed — the WCC test had dead code (built a sole-`y` surrogate but never
  asserted); replaced with real matched-null tests on `.mv_wcc_cell`/`.mv_wdtw_cell`/
  `.mv_granger_cell` (sole surrogate = `y` ⇒ `null_mean == observed`, both Granger directions),
  polarity tests retained. (2) AC8 validation strengthened — `sim_dyad` regression now pins the
  recommendation to the recovered 0.5 Hz cycle (`window_sec == 2`, `window_size == 160`,
  `sig_rate == 1`, `median_es > 1.5`); a true uncoupled negative control (i.i.d. noise dyads) was
  added asserting the detectability gate fails (warning fired, `sig_rate < 0.5`). (3) `n_cells`
  naming disambiguated — `robustness` now carries `n_cells` (grid total) **and** `n_valid`
  (computable cells); `print`/plot title/`glance()` updated (`glance` gains an `n_valid` column;
  `pct_significant` is over `n_valid`). (4) `wcc-params` vignette dropped (superseded by
  `choosing-parameters`; overlap discussion folded in, `wcc-workflow` link redirected). (5) clean
  code-line wraps reduced line-length lints; remaining lints are cli strings/roxygen/idiomatic
  predicates (package-wide style) and glue-variable `object_usage` false positives.

## Current focus

**Post-release maintenance and next-estimator work.** M7 (`v0.1.0`) is complete and ready to submit
to CRAN (human-gated — see `cran-comments.md` pre-submission checklist). The next milestone is M8.

- **M8 — Phase synchrony estimator (next).** Add a windowed phase-synchrony estimator to the
  `bsync_surface` framework established in M5. See `DESIGN.md` §15 for scope.

See `DESIGN.md` §15 for the full roadmap (M8 phase synchrony; M9 wavelet coherence; M10 CRQA/MEA;
M11 group-level workflow; M12 expanded surrogates).

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
`dplyr`, `future.apply`, `generics`, `ggplot2`, `grDevices`, `gsignal`, `rlang`, `scales`,
`tibble`, `utils`. `Suggests`: `future`, `knitr`, `rmarkdown`, `spelling`, `testthat (>= 3.0.0)`,
`vdiffr`. `generics` and `tibble` added in M5 for the tidy interface. **Do not grow `Imports`
without flagging it** — prefer base R or the existing stack, and put heavy compute in the C++ cores,
not new R deps.

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
- **tidy/glance/as_tibble methods** — built in **M5** (done).
- **IAAFT / segment-shuffling / pseudo-dyad surrogates** — resolved to add (DESIGN.md §6/§14),
  scheduled **M12**; the pseudo-dyad (between-dyad rMEA) generator depends on M11's `dyad_list`.
- **A unified `bsync_ts` preprocessing object** — logged in DESIGN.md §14; specify before building.
- **OpenMP** — not adopted until the M2 decision; no `#pragma omp` ships before then.

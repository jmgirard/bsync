# M014: IAAFT surrogate generator

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Resolves:** —
- **Surface tier:** user-facing — adds an exported generator and a new `surrogate_method` value
- **Branch/PR:** m014-iaaft-surrogates

## Goal

Add an exported IAAFT surrogate generator that keeps the value distribution and the power spectrum of a series. Offer it as a `surrogate_method` in `synchrony_multiverse()` and `autotune_wcc()`.

## Scope

**In:**
- `generate_surrogate_iaaft(y, n_surrogates = 100, max_iter = 1000)` in `R/surrogate_generation.R`, after schreiber1996 (p. 2). Each column starts from a random shuffle of `y` without replacement. Each iteration replaces the Fourier amplitudes with those of `y` and keeps the phases. It then transforms back and rank-orders the result, so that the result takes exactly the values of `y`. When the rank-ordering no longer changes a column, that column stops (p. 2). At `max_iter`, every column stops. `match = "spectrum"` (the default) returns the last spectrum-adjusted series of each column: its Fourier amplitudes equal those of `y`, and its values are close to the values of `y`. `match = "values"` returns the rank-ordered series, which holds exactly the values of `y`.
- A warning that names how many columns reached `max_iter` without converging. The paper advises reporting a failure when the iteration does not reach its target (p. 2). The matrix is still returned.
- An error for `NA` or non-finite values in `y`, with a message that points to `impute_ts_gaps()`. An error when `length(y) < 3`. Errors for invalid `y`, `n_surrogates`, or `max_iter`.
- `"iaaft"` as a `surrogate_method` value in `synchrony_multiverse()` and `autotune_wcc()`. Like `"phase"`, it is generated once per call and reused across cells.
- Roxygen with the null this tests, the `max_iter` rationale, and the citation. A `_pkgdown.yml` row, an IAAFT section in the surrogate-testing vignette, NEWS, `inst/WORDLIST`.
- Source note `cairn/references/schreiber1996.md` and DESIGN §2, §6, §13, §14 #8 updates.

**Out:**
- Segment-shuffling surrogates → M015.
- The AAFT algorithm and a windowed or smoothed target spectrum (schreiber1996, p. 3 remarks): not asked for, no row.
- A C++ core for the iteration. The work log names the measurement that justifies one.
- Pseudo-dyad surrogates in the multiverse: their existing candidate row.

## Acceptance criteria

- [ ] AC1: Tests on seeded inputs show that `generate_surrogate_iaaft()` returns a numeric matrix with `length(y)` rows and `n_surrogates` columns for both `match` values. With `match = "values"`, `sort(surr[, j])` is identical to `sort(y)` for every column `j`. With `match = "spectrum"`, `Mod(fft(surr[, j]))` equals `Mod(fft(y))` (`expect_equal`) for every column `j`. The inputs include an even length, an odd length, and a series with tied values.
- [ ] AC2: A test file holds a plain reimplementation of the schreiber1996 (p. 2) iteration. It uses explicit loops and the DFT as the explicit sum in the paper's S_k formula. It uses the same stopping rule and the same sequence of `sample()` draws. The test uses at least two seeded series of length 32 or less, one even and one odd. For each, every column of `generate_surrogate_iaaft()` equals the reference output for both `match` values (`expect_equal`).
- [ ] AC3: A test uses a seeded series where no column reaches `max_iter`. With `match = "values"`, it shows that one more iteration of the reference step leaves every output column unchanged. It also computes the spectral discrepancy of schreiber1996 (p. 3), unsmoothed. For each column, this discrepancy is smaller than that of a random shuffle of `y` drawn in the test. With `match = "spectrum"` under the same seed, rank-ordering each column with `sort(y)[rank(col, ties.method = "first")]` gives the `"values"` column exactly. On a seeded AR(1)-cube series of length 128, each `"spectrum"` column has a smaller `mean(abs(sort(col) - sort(y)))` than the mean of that distance over 20 phase surrogates of `y` drawn in the test.
- [ ] AC4: A seeded simulation test (`skip_on_cran()`) draws at least 400 independent pairs of length 512. Each partner is an AR(1) series x_n = 0.7 x_{n-1} + eta_n, observed as s_n = x_n^3 (schreiber1996, p. 2). The test runs `wcc_surrogate()` with `window_size = 32`, `lag_max = 4`, `window_increment = 16`, and at least 19 IAAFT surrogates per pair at the default `match`. The proportion of pairs with p <= .05 is at most 0.08.
- [ ] AC5: A test calls the generator with `max_iter = 1` on a series that does not converge in one iteration, once for each `match` value. Each call warns with a message that names the number of unconverged columns, and it returns the full matrix. Tests fire each abort branch and assert each message with a regexp. The branches are `NA` in `y`, a non-finite value, non-numeric `y`, and `length(y) < 3`, with a test that length 3 runs. The others are an invalid `n_surrogates`, an invalid `max_iter`, and an invalid `match`.
- [ ] AC6: A test runs a one-cell `synchrony_multiverse(estimator = "wcc", surrogate_method = "iaaft")`. Its `p` equals the `p` of `wcc_surrogate()` on `generate_surrogate_iaaft()` output drawn after the same `set.seed()`, with the cell's window, lag, and increment in samples. A test shows that a `c("phase", "iaaft")` grid has both values in its `surrogate_method` column. A test shows that `autotune_wcc(surrogate_method = "iaaft")` returns a result whose grid holds `"iaaft"`.
- [ ] AC7: The roxygen of `generate_surrogate_iaaft()` states the null it tests and the default `match`. For each `match` value, it states which property is exact and which is close. It states the `max_iter` default and its reason, and the `NA` policy. It cites Schreiber & Schmitz (1996) with its DOI. It has a runnable `sim_dyad` example that shows the exact property of the default form. The `surrogate_method` docs of both functions name `"iaaft"`. The surrogate-testing vignette has an IAAFT section. NEWS.md has an entry. `pkgdown::check_pkgdown()` passes. `devtools::document()` leaves no uncommitted diff in `man/` or `NAMESPACE`. `devtools::test()` passes. `devtools::check()` gives 0 errors, 0 warnings, and no NOTE that `main` does not also give.

## Coverage

- AC1 → T2, T3
- AC2 → T2, T3
- AC3 → T2, T3
- AC4 → T4
- AC5 → T2, T3
- AC6 → T5
- AC7 → T1, T3, T6, T7

## Tasks

- [x] T1: Write `cairn/references/schreiber1996.md` from the source-note template. The source is the arXiv copy on the shelf (`sources/schreiber1996.pdf`, chao-dyn/9909041, PRL 77, 635). Extract the iteration and stopping rule (p. 2), the discrepancy measure and the AR(1) cube process (pp. 2-3), and the convergence remarks (p. 3). Add the INDEX.md line.
- [x] T2: Write the AC1, AC2, AC3, and AC5 tests first. Put the plain reference in the test file with an oracle provenance header (DESIGN §13). Then implement `generate_surrogate_iaaft()` in R next to `generate_surrogate_phase()` (`R/surrogate_generation.R:84`). Use `stats::mvfft` over the columns that have not converged. Time it on `sim_dyad` (2400 samples, 100 surrogates) and record the time in the work log.
- [x] T3: Add `match = c("spectrum", "values")` to `generate_surrogate_iaaft()` and to the plain reference, and extend the AC1, AC2, AC3, and AC5 tests to both forms. Make the warning hint name the property that is close for the given `match`, and assert it in both AC5 calls. Rewrite the source note's bsync-choices bullet to match.
- [x] T4: Write the AC4 calibration test with the design AC4 names. Record the observed rejection rate and the run time in the work log.
- [x] T5: Add `"iaaft"` to the `surrogate_method` choices and the generation branch in `synchrony_multiverse()` (`R/multiverse.R:192`, `R/multiverse.R:279`). Add it to the `autotune_wcc()` docs (`R/autotune.R:80`). Write the AC6 tests.
- [x] T6: Write roxygen with `@md`, `@references`, and `@seealso` links both ways with the other generators. The roxygen gives the reason for the default, the measured `"values"` rejection rate, and "close" as measured (closer to the values of `y` than a phase surrogate, on a skewed AR series). Reword DESIGN §6 so that each `match` form keeps one property exactly and the other closely. Add the `_pkgdown.yml` row, the vignette section, the NEWS entry, and the WORDLIST words. Update the method count in section 1 of the vignette. Update DESIGN §2, §6, §13 (oracle records), and §14 #8. Run `devtools::document()`.
- [x] T7: Run `devtools::test()`, `spelling::spell_check_package()`, `pkgdown::check_pkgdown()`, and `devtools::check()`. Run `air format` only on files that were air-clean before the edit.

## Work log

- 2026-10-06: created by /milestone-plan. It takes the IAAFT half of the "Expanded surrogate generators" candidate row. M015 takes the segment-shuffling half.
- 2026-10-06: question set: can the agent fetch the primary sources itself — yes. The arXiv copy of schreiber1996 is on the shelf. Nothing was installed.
- 2026-10-06: question set: offer the new generators in `synchrony_multiverse()` and `autotune_wcc()` — yes, both.
- 2026-10-06: plan chose R over a C++ core, because `stats::mvfft` transforms all columns at once and the rank step is a native sort; falsified by a timing where IAAFT generation takes longer than the cell analysis of a typical multiverse run.
- 2026-10-06: plan chose the fixed-point stop over a spectral tolerance, because the paper says the iteration ends at a fixed point (p. 2) and it needs no new argument; falsified by typical series that do not converge within 1000 iterations.
- 2026-10-06: plan chose to return the rank-ordered series over the last spectrum-matched series, after the paper (p. 2); falsified by a user need for an exact spectrum, which `generate_surrogate_phase()` already gives.
- 2026-10-06: plan chose to abort on `NA` over imputing inside the generator, because the FFT cannot take `NA` and DESIGN forbids silent changes; falsified by none expected.
- 2026-10-06: oracle plan: closed-form (plain reimplementation, AC2), invariant (AC1, AC3, AC6), and simulation-coverage (AC4). DESIGN §13 layer 3 validates generators by properties, so there is no frozen pin.
- 2026-10-06: criteria audit (full mode, fresh Opus reader): 2 findings, both fixed. AC5 now names the length threshold (`length(y) < 3`). AC7 no longer binds the Review record. The audit also corrected the Scope reading of the paper: the paper advises reporting a failure and names a spectral target as its main stop, and the fixed point is where the iteration ends. AC1 to AC4 and AC6 passed.
- 2026-10-06: T1 done: wrote `cairn/references/schreiber1996.md` and its INDEX line from the arXiv copy, pp. 1-4.
- 2026-10-06: T2 done: `generate_surrogate_iaaft()` and `tests/testthat/test-surrogate-iaaft.R`. 100 surrogates of `sim_dyad$x_B` (2400 samples) took 1.0 s elapsed, with no unconverged column. Planted defects: squared target amplitudes and a reversed shuffle both turn the reference test red. A changed tie rule in the rank step does not, because the adjusted series has no ties on continuous data. Full suite: 271 tests, 0 failed.
- 2026-10-06: T2 choice: the reference test uses the paper's +i DFT sign while the generator uses R's -i sign, so the match also shows the result does not depend on the sign convention.
- 2026-10-06: T3 calibration (old T3) found the rank-ordered output too liberal: 400 null pairs, 19 surrogates, window 32, lag 4, increment 16 gave 10.0% (n = 128) and 10.75% (n = 512) at p <= .05. Lag-1 autocorrelation was 0.698 in the data and 0.635 in the surrogates. The last spectrum-adjusted series gave 6.5% and 5.25%. Phase gave 6.0% and 1.5%, circular 7.0% and 4.25%.
- 2026-10-06: stop at the user: which series to return. The user chose to add `match = c("spectrum", "values")` with the spectrum form as the default.
- 2026-10-06: substantive amendment: Scope bullet 1 last sentence, AC1, AC2, AC3, AC4, AC5, AC7 first part reworded for `match`. A new T3 holds the `match` work. The old T3 to T6 became T4 to T7, and Coverage moved with them. The Goal wording stays: each form keeps one property exactly and the other closely.
- 2026-10-06: re-audit: AC1 (full) — the spectrum form's values-closeness had no criterion; fixed by the AC3 addition.
- 2026-10-06: re-audit: AC2 (full) — nothing.
- 2026-10-06: re-audit: AC3 (full) — no check that both forms come from one iteration, and no closeness check; both added.
- 2026-10-06: re-audit: AC4 (full) — 200 pairs and 0.09 did not separate the forms; fixed to 400 pairs of length 512, the named design, and 0.08.
- 2026-10-06: re-audit: AC5 (full) — the warning hint was wrong for the spectrum form; fixed by asserting the warning for both forms.
- 2026-10-06: re-audit: AC7 (full) — the roxygen did not need to say which property is exact per form; fixed.
- 2026-10-06: re-audit: AC1 (full) — nothing.
- 2026-10-06: re-audit: AC2 (full) — nothing.
- 2026-10-06: re-audit: AC3 (full) — "close" is checked on one series type; handled in T6 by stating "close" as measured, wording unchanged.
- 2026-10-06: re-audit: AC4 (full) — nothing. P(more than 32 of 400) is 0.0076 at 5.25% and 0.897 at 10%.
- 2026-10-06: re-audit: AC5 (full) — the hint text was not asserted; handled in T3 by asserting the hint for both forms, wording unchanged.
- 2026-10-06: re-audit: AC7 (full) — the roxygen did not need to give the default's reason and the `"values"` rate; handled in T6, wording unchanged. Also: the source note and DESIGN §6 still state the old choice; handled in T3 and T6.
- 2026-10-06: T3 done: `match = c("spectrum", "values")` via `rlang::arg_match()`, the reference returns both forms, tests cover both, and the warning hint names the close property per form. The roxygen was rewritten in T3 because it stated the old default; T6 adds the measured rates. Full suite: 273 tests, 0 failed. The two examples print `TRUE` for the exact property of each form.
- 2026-10-06: T4 done: `tests/testthat/test-surrogate-calibration.R`, seed 20261006. The spectrum form rejected 17 of 400 pairs (4.25%) at p <= .05. The values form rejected 21 of the first 200 pairs (10.5%), against 4.5% for the spectrum form on the same 200 pairs; the test asserts both. The file runs in 51 s. Full suite: 274 tests, 0 failed.
- 2026-10-06: T4 speed change: the rank step now ranks all active columns in one stable `order()` call in place of `apply(rank)`. Output was identical (`identical()`) on four inputs for both forms, and 20 calls at n = 512 went from 1.92 s to 1.22 s.
- 2026-10-06: T5 done: `"iaaft"` in the `synchrony_multiverse()` choices and generation branch, and in both functions' docs. Neither roxygen block uses markdown, so the new text uses Rd markup (`\code{\link{}}`) for the link to render. Three tests in `test-multiverse.R` cover AC6. Full suite: 277 tests, 0 failed.
- 2026-10-06: T6 done: roxygen gives the measured rates with their design and date, the vignette has section 1.3 IAAFT (pseudo-dyads became 1.4), plus the `_pkgdown.yml` row, NEWS, the README reference, three WORDLIST words, and DESIGN §2, §6, §13, §14 #8. `spelling::spell_check_package()` finds 0 words. `pkgdown::check_pkgdown()` finds no problems. `build_readme()` also rewrote `man/figures/README-example-plot-1.png`, which this change does not affect, so that file was restored.
- 2026-10-06: T7 in progress: `air format` left the new and air-clean files unchanged. Two lintr line-length hits of this branch were fixed; the other hits are in older code. `devtools::check()` is running.
- 2026-10-06: claim audit: 38 claims read, 2 corrected — `R/surrogate_generation.R` (the Fig. 2 page cite is p. 2, the 1/i statement p. 3), `tests/testthat/test-surrogate-iaaft.R` (`spectral_discrepancy()` now uses the paper's amplitude scale; the AC3 test still passes). The reader's re-read of the two corrections is pending.
- 2026-10-06: claim audit re-read: both corrections OK. The reader also flagged one 82-character roxygen line, rewrapped in `7f31563`; the IAAFT block now has no lint hits.
- 2026-10-06: T7 done: `devtools::check()` on `7f31563` gave 0 errors, 0 warnings, 0 notes (2m 14s). Status set to review.

## Decisions

- 2026-10-06: `generate_surrogate_iaaft()` returns the spectrum-adjusted series by default (`match = "spectrum"`), with the rank-ordered series as `match = "values"`. This reverses the plan's choice of the rank-ordered series. In a WCC size simulation, the rank-ordered form rejected 10.0% and 10.75% of true-null pairs at p <= .05, and the spectrum form rejected 6.5% and 5.25% (work log above). The rank step lowers the autocorrelation, so a statistic that grows with autocorrelation reads high against those surrogates. The user chose this at the implement stop.

## Review

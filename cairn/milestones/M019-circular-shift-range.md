# M019: Circular surrogates in the multiverse: a shift range per cell

- **Status:** review
- **Priority:** high
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Resolves:** —
- **Surface tier:** user-facing — changes the results and warnings of the exported `synchrony_multiverse()` and `autotune_wcc()`
- **Branch/PR:** m019-circular-shift-range

## Goal

Give each circular cell of `synchrony_multiverse()` its own range of valid circular shifts, with distinct shifts only, and skip a cell that has too few shifts to reach p <= .05.

## Scope

**In:** In `synchrony_multiverse()`, a wgranger circular cell shifts `y` by at least one window (D-005). A wcc or wdtw circular cell shifts `y` by at least `2 * lag_max` for its own lag, with one surrogate matrix per distinct lag. Each cell uses `min(n_surrogates, k)` distinct shifts, where k is its number of valid shifts. A cell that fits its window and has k < 19 is skipped with a warning. This rule replaces the M018 rule (no valid shift at `4 * lag_max`). The docs, messages, and NEWS of `synchrony_multiverse()` and `autotune_wcc()` change to match. CLAUDE.md invariants 2, 6, and 7 hold: the null statistic does not change, nothing reseeds, and the result stores no shifts.

**Out:** `generate_surrogate_circular()` keeps its output and its sampling with replacement (no change planned). The wphase segment and IAAFT rows stay candidate rows. Pseudo-dyad cells in the multiverse stay a candidate row.

## Acceptance criteria

- [x] AC1: For `estimator = "wgranger"`, `synchrony_multiverse()` draws the surrogates of each circular cell from circular shifts s with `window_size <= s <= n - window_size`. Here n is `length(y)`, and `window_size` is the window of the cell in samples. A test captures the surrogate matrix that each computed Granger circular cell receives and recovers the shift of each column by matching it to the rotations of `y`. The test uses `window_sec = c(2, 5)` at `sample_rate = 10` on series of length 200, 201, 202, 203, and 204. It asserts that range for each computed cell, and that none of the five calls stops with an error.
- [x] AC2: For `estimator = "wcc"` and `"wdtw"`, `synchrony_multiverse()` draws the surrogates of each circular cell from circular shifts s with `2 * lag_max <= s <= n - 2 * lag_max`. Here `lag_max` is the lag of that cell in samples, not the largest lag of the grid. The test series has length 401, with `sample_rate = 1`, `window_sec = c(10, 200)`, `lag_sec = c(5, 20, 100)`, and `n_surrogates = 50`. For wcc and for wdtw, the test asserts these results. The cell with window 10 and lag 5 uses 50 distinct shifts in `10:391`, and at least one of them is outside `200:201`. The cell with window 200 and lag 20 uses 50 distinct shifts in `40:361`. The cell with window 200 and lag 100 has k = 2 and is skipped with a `bsync_circular_skipped` warning.
- [x] AC3: A computed circular cell, for any estimator, uses `min(n_surrogates, k)` surrogates. Here k is the number of valid shifts of the cell under AC1 or AC2. No two surrogates of one cell use the same shift. The `n_surrogates` column of the grid holds that count. For wcc, wdtw, and wgranger, a test asserts the count and the distinct shifts of a cell with 19 <= k < `n_surrogates`.
- [x] AC4: A circular cell that fits its window and has k < 19 is skipped, because no add-one p-value of its surrogates can reach p <= .05. A cell fits when n >= `window_size + 2 * lag_max`, or n >= `window_size` for wgranger. The skip also applies when the call asks for fewer than 19 surrogates. It gives a `bsync_circular_skipped` warning that names the cell and states k. The grid row of the cell has `NA` results and `n_surrogates` 0. A cell too short for its window keeps the too-short handling and gives no `bsync_circular_skipped` warning. For wcc and for wgranger, tests assert that a cell with k = 18 is skipped, that a cell with k = 19 is computed, and that a cell too short for its window gives no circular warning.
- [x] AC5: The pinned call is `synchrony_multiverse(sim_dyad$x_A, sim_dyad$x_B, estimator = "wcc", sample_rate = 80, window_sec = c(1, 2), lag_sec = 0.5, surrogate_method = c("phase", "circular", "segment"), n_surrogates = 19)` after `set.seed(20261011)`. Its circular cells share one lag, and each has k >= 19. For every cell of that call, `observed`, `null_mean`, `null_sd`, `es`, `p`, and `n_windows` equal main at c13ffca with difference 0. `generate_surrogate_circular()` returns the same matrix and warnings as at c13ffca for the same seed in three pinned cases: `lag_max = NULL`, `lag_max = 5` on a long series, and a `lag_max` that leaves fewer valid shifts than `n_surrogates`.
- [x] AC6: The `surrogate_method` docs of `synchrony_multiverse()` and `autotune_wcc()` state the shift ranges of AC1 and AC2, the count of AC3, and the skip of AC4, including the skip below 19 requested surrogates. The `bsync_circular_skipped` messages of both functions state the k < 19 rule. A sweep reads every hit of `grep -n -C3 -e circular -e '4 \* lag_max' -e 'n / 4'` over `R/*.R` and `vignettes/*.Rmd`. The sweep finds no doc or message text that states an old rule: one shared circular matrix, the largest lag of the grid, `n / 4` for wgranger, or a skip only at `4 * lag_max`. NEWS.md states which results change. These are circular cells of wgranger grids, of wcc and wdtw grids with more than one lag, and of cells with fewer valid shifts than `n_surrogates`. Other cells of a call with circular cells can also change, because the random draws shift. NEWS also states the new skip.

## Coverage

- AC1 → T2, T3
- AC2 → T2, T3
- AC3 → T2, T3
- AC4 → T2, T3, T4
- AC5 → T1, T4
- AC6 → T5

## Tasks

- [x] T1: In `R/surrogate_generation.R`, add an internal helper that builds circular surrogates from a minimum shift and a number of draws. Make `generate_surrogate_circular()` call it with `2 * lag_max` and keep its sampling with replacement and its warning. First, derive the three AC5 generator pins from c13ffca with `git archive` into the scratchpad (LESSONS, M018). Then add the pin tests.
- [x] T2: In `R/multiverse.R` (about lines 423-485), give each circular cell its own shift range. Build one matrix per distinct lag (wcc, wdtw) or per distinct window (wgranger), each with `min(n_surrogates, k)` distinct shifts. Keep these draws in the circular slot of the generation order, so a single-lag call draws its random numbers as main does (AC5). Skip cells that fit and have k < 19, and state k in the warning. Remove the M018 shared-lag logic. Set `n_surrogates` per cell.
- [x] T3: Add the tests of AC1 to AC4. Mock `.mv_wcc_cell`, `.mv_wdtw_cell`, and `.mv_granger_cell` to capture `y_surr`, and recover each shift by matching the columns to rotations of `y`.
- [x] T4: Update the tests that relied on circular cells with few shifts. These include `test-multiverse.R` near lines 255, 605, 623, 791, 901, 926, and 953, and `test-autotune.R` near line 473. That autotune test needs a series where the lag-2 cell has k >= 19, so it still selects lag 2. Add the AC5 multiverse pin. Re-derive every pin in `test-multiverse-pins.R` from a call with circular cells, segment pins included, and record the reason in the file header. Make sure that the vdiffr snapshots do not change (their fixtures use phase).
- [x] T5: Update the roxygen and cli messages in `R/multiverse.R` and `R/autotune.R` (about lines 85-95 and 305-320), and add the NEWS entry. Run the AC6 sweep. Run `devtools::document()`, `spelling::spell_check_package()`, `air format` on air-clean files, and `devtools::check()`.
- [x] T6: Review return 1. A cell with lag 0 (a 1-sample window) gets minimum shift 0 and stops the call. Floor the wcc and wdtw minimum shift at 1 sample, state it in the roxygen, and add a regression test that fails before the fix.

## Work log

- 2026-10-11: created by /milestone-plan. Absorbs the candidate rows "Granger circular surrogates in the multiverse" and "Shared circular matrix with few shifts". Inbox sweep: 0 open issues, 0 open PRs.
- 2026-10-11: question set: scope is circular surrogates in the multiverse. Granger rule is at least one window. For wcc and wdtw, one matrix per lag. Cells with fewer than 19 shifts are skipped with a warning. The user approved these results changes at the gate.
- 2026-10-11: plan gate chose a minimum Granger shift of one window over 2 * ar_order, because a shift of 2 or 3 samples keeps much of the coupling in smooth signals. Falsified by a power simulation where the window rule loses power against the ar_order rule at equal size.
- 2026-10-11: plan gate chose one circular matrix per lag over one shared matrix for the largest lag, because the shared matrix leaves small-lag cells with few shifts. Falsified by evidence that cells of one grid must share their surrogate draws to be compared.
- 2026-10-11: plan gate chose to skip cells with k < 19 over computing them with a warning, because it matches the segment rule (D-002, D-004). Falsified by users who need such cells counted in the robustness summaries.
- 2026-10-11: criteria audit (full mode, fresh Opus reader) returned 9 findings, all fixed. The AC2 test now separates new code from main. Tests read shifts by capturing `y_surr`. AC1 probes two windows. D-005 records the Granger rule. AC5 is narrowed to the pinned call and three generator cases. NEWS covers changed non-circular cells. The AC6 sweep reads hits with context and covers autotune messages. AC4 keeps too-short cells apart and states the skip below 19 requested surrogates. T4 names the autotune test outcome. AC3 probes wdtw.
- 2026-10-11: lessons applied: M018 `git archive` pin derivation (T1, T4), M017 snapshots (T4), M010 air on clean files only (T5), M011 spelling (T5).
- 2026-10-11: implement started on branch m019-circular-shift-range.
- 2026-10-11: T1 done. Added internal `circular_shift_matrix()` and `draw_circular_shifts()` in `R/surrogate_generation.R`. `generate_surrogate_circular()` builds its matrix with the first. Three pins from c13ffca (seeds 101-103) pass in `test-surrogate.R`. `draw_circular_shifts()` uses `x[sample.int()]`, which is what `sample()` does for length > 1, so it draws as main does when k >= n_surrogates.
- 2026-10-11: T1 full suite: 2678 expectations, 0 failed.
- 2026-10-11: T2 done. `synchrony_multiverse()` computes `circ_min_shift` (2 * lag_max, or window for wgranger), `circ_k`, and `circ_skip` per cell, and draws one matrix per distinct minimum shift in the circular slot of the method order. The skip warning carries `window_size`, `lag_max`, and `n_shifts` fields. Implementation choice: grid order of first occurrence for the matrices (a descending order keeps one more multi-lag cell equal to main, but no criterion asks for it).
- 2026-10-11: T3 done. New `test-multiverse-circular.R` (7 tests) recovers each shift from the captured `y_surr`. Against c13ffca, tests 1 to 6 fail and test 7 (too-short cells give no circular warning) passes, as it guards unchanged behavior.
- 2026-10-11: T4 done. Moved the two one-window tests in `test-multiverse.R` to phase surrogates, since their point is the window count. Rebuilt the short-series circular test on 30 samples. Mocked `circular_shift_matrix()` in the two-matrices test. The autotune circular test uses 30 samples, so lag 2 (23 shifts) runs and lag 5 (11 shifts) is skipped. Re-derived cells 4 to 9 of the two-lag pin (circular and segment changed, phase did not). The single-lag pin uses hexadecimal doubles, because R parsed the decimal literal 0.00076342284885224038 one ulp below the value that main and the branch both compute. It passes on the branch and on c13ffca. No vdiffr snapshot changed. Full suite: 3266 expectations, 0 failed.
- 2026-10-11: T5 done. Rewrote the circular text of the `surrogate_method` roxygen of `synchrony_multiverse()` and `autotune_wcc()`, and the autotune skip message. Replaced the M018 NEWS bullet on circular skips with four bullets (shift range, distinct shifts, skip below 19, which results change). The AC6 sweep read every hit with context in `R/*.R` and `vignettes/*.Rmd`. One stale code comment (`R/multiverse.R` near line 647) was fixed. No vignette states an old rule. `document()` rewrote two Rd files, spelling found no errors, air ran on the three air-clean files, and `check()` gave 0 errors, 0 warnings, 0 notes.
- 2026-10-11: claim audit: 41 claims read, 4 corrected — NEWS.md, R/multiverse.R (fixed in d37c5ed; the same reader re-read all 4 as correct). Reader note left as is: when every circular cell of a dyad is too short for its window, main still drew a circular matrix and the branch does not, so a later dyad of `autotune_wcc()` can draw different numbers.
- 2026-10-11: after the corrections, `document()` rewrote one Rd, the full suite gave 3266 expectations with 0 failed, and spelling found no errors. Status set to review.
- 2026-10-11: review return 1: diff-bug #1. A window that rounds to 1 sample caps its lag to 0, so its circular cell gets minimum shift 0 and draws shift 0 or n. `circular_shift_matrix()` then stops the whole `synchrony_multiverse()` call ("number of items to replace is not a multiple of replacement length"; n = 100, `window_sec = c(1, 10)`, `lag_sec = 2`, wcc and wdtw). Main returned results for the same call. The k formula also counts n + 1 shifts for m = 0. Fix the minimum shift for lag 0, and add a regression test that fails before the fix.
- 2026-10-11: review fix-now items committed before the return: `n_surrogates` docs and comment (f8bd33a); `glance()` docs, the M018 NEWS bullet, the zero-shift skip message, and the NEWS note on shifted draws (review-fix commit). Status set to in-progress.
- 2026-10-11: implement resumed for review return 1 and added T6 (minor amendment).
- 2026-10-11: T6 done. The wcc and wdtw minimum shift is now `pmax(1L, 2L * l_samp_all)`, and the roxygen says "and at least 1". A new test in `test-multiverse-circular.R` (n = 100, `window_sec = c(1, 10)`, `lag_sec = 2`, wcc and wdtw) failed before the fix with "number of items to replace is not a multiple of replacement length". After the fix it passes: the lag-0 cell uses shifts 1 to 99, and the lag-2 cell computes with 93 shifts. AC2 still holds as written, because shifts of at least 1 meet `2 * lag_max <= s` at lag 0. Full suite: 3476 expectations, 0 failed.
- 2026-10-11: claim audit: 14 claims read, 2 corrected — NEWS.md, tests/testthat/test-multiverse-circular.R (lines added since d37c5ed; one optional comment clarification in R/multiverse.R also taken; fixed in d4b0a88). On the one re-read, the NEWS sentence on main's old draws was still wrong (main drew for a fitting smaller-lag cell, and wgranger could stop with an error). With no second pass allowed, the sentence was narrowed to the verified part: a call in which every circular cell is skipped draws no circular surrogates, so its later draws can differ.
- 2026-10-11: spelling clean, `document()` no change. Status set to review.

## Decisions

- 2026-10-11 (review return 1): the wcc and wdtw minimum circular shift is `max(1, 2 * lag_max)`. A 1-sample window caps its lag to 0, and shifts of 0 or n leave `y` in place. Shifts of at least 1 still meet D-005's bound `2 * lag_max <= s`, so D-005 stands unchanged.

## Review

Review on 2026-10-11, branch m019-circular-shift-range at b6001ec, main e3d98b9 unmoved.

- AC1: `test-multiverse-circular.R` test 1 (230 expectations, 0 failed) runs wgranger with `window_sec = c(2, 5)` at 10 Hz on n = 200 to 204. It recovers every shift of both windows from the captured `y_surr`, asserts `window_size <= s <= n - window_size` for each, and asserts no error. The same test failed on c13ffca during implement.
- AC2: test 2 (226 expectations, 0 failed) runs wcc and wdtw on n = 401 with windows 10 and 200, lags 5, 20, and 100, and 50 surrogates. Cell (10, 5) uses 50 distinct shifts in `10:391`, with at least one outside `200:201`. Cell (200, 20) uses 50 distinct shifts in `40:361`. Cell (200, 100) has k = 2 and gives one `bsync_circular_skipped` warning with `n_shifts` 2, `NA` es, and `n_surrogates` 0.
- AC3: test 3 (99 expectations, 0 failed) gives k = 30 below 50 requested surrogates. For wcc and wdtw (n = 49, window 10, lag 5), the sorted shifts equal `10:39`. For wgranger (n = 69, window 20), they equal `20:49`. Both are 30 distinct shifts, and the grid `n_surrogates` is 30 in each case.
- AC4: tests 4 to 7 (28 expectations, 0 failed). For wcc (window 4, lag 1), n = 21 gives k = 18: skipped with one warning, `n_shifts` 18, message "18 distinct shifts", `NA` es and p, `n_surrogates` 0. n = 22 gives k = 19 and computes with 19. For wgranger (window 10), n = 37 is skipped with k = 18, and n = 38 computes. With 9 surrogates requested, k = 18 is still skipped, and k = 27 computes with 9. A cell too short for its window gives the too-short message and no circular warning, for wcc and for wgranger. A direct run of wgranger with n = 37 printed a warning that names the cell ("Skipped a circular cell with `window_size` = 10 and `window_increment` = 1 samples") and states "18 distinct shifts".
- AC5: `test-multiverse-pins.R` test 3 (5 expectations, 0 failed, tolerance 0) passes on the branch and, rerun now, on c13ffca. So every cell of the pinned single-lag call equals main. `test-surrogate.R` tests 8 to 10 pass. A direct run of the three generator cases (seeds 101 to 103) on the branch and on c13ffca gave `identical()` TRUE for the matrices and the warnings together. The third case gives the replacement warning on both.
- AC6: `man/synchrony_multiverse.Rd` (`surrogate_method`) states m = 2 * lag_max or window_size, the range m to n - m, k = n - 2 * m + 1, `min(n_surrogates, k)` distinct shifts, and the skip below 19, also below 19 requested surrogates. `man/autotune_wcc.Rd` states the wcc form and the same skip. Both skip messages state the 19-shift rule. The rerun sweep read the context of every hit and found 3 hits on the old-rule words, none stating an old rule. Two state the new rule (`R/autotune.R:91`, `R/multiverse.R:203`), and one is the standalone generator advice in `vignettes/wphase-workflow.Rmd:146`. NEWS.md lists the changed results (wgranger, multi-lag wcc and wdtw, k below `n_surrogates`, and later-drawn cells of other methods) and the new skip.

Consistency gate (pass 1): `cairn_validate` all checks passed. `document()` made no diff. README not touched. `pkgdown::check_pkgdown()` found no problems. NEWS has the entry. No new top-level files. `check()` gave 0 errors, 0 warnings, 0 notes.

spawned: diff-bug, blame-history, prior-review

- diff-bug #1: a 1-sample window caps its lag to 0, so its circular cell gets minimum shift 0, draws shift 0 or n, and `circular_shift_matrix()` stops the whole call (n = 100, `window_sec = c(1, 10)`, `lag_sec = 2`, wcc and wdtw). Main returned results — floor return (load-bearing crash), review return 1.
- diff-bug #2: the `n_surrogates` grid-column docs still say circular surrogates repeat shifts — fix now, fixed f8bd33a (same as prior-review #1).
- diff-bug #3: the `glance()` docs and the M018 NEWS bullet name only segment cells as using fewer surrogates — fix now, fixed in the review-fix commit.
- diff-bug #4: the aggregated `bsync_circular_skipped` warning of `autotune_wcc()` has none of the fields that the `synchrony_multiverse()` warning carries — follow-up, new candidate row "Circular-skip condition fields".
- diff-bug #5: with k = 0 the skip message stated inverted bounds ("at least 30 and at most 10") — fix now, the message now says no valid shift exists. Fixed in the review-fix commit.
- diff-bug #6: the wgranger skip warning sets `lag_max` to `NA` with no docs — follow-up, same row as #4.
- diff-bug #7: a call where no circular cell is computed no longer draws circular surrogates, so later draws shift, and NEWS did not say so — fix now, NEWS sentence added in the review-fix commit.
- diff-bug #8: the re-derived two-lag pins guard only later drift — reject, planned change (T4 re-derives them, and AC2 tests the ranges).
- diff-bug #9: no test covers a lag-0 or 1-sample circular cell — goes with review return 1 (its regression test).
- diff-bug #10: no test asserts two matrix builds for two lags, or independent draws across estimators — reject, false: the AC2 test recovers each cell's shifts and fails if two lags share a matrix.
- diff-bug #11: the draw equivalence rests on `sample()` calling `sample.int()` — reject, false as a defect: the comment states it, and the AC5 pins fail if it stops holding.
- blame-history #1: cells with 1 to 18 shifts that M018 computed are now skipped — reject, planned change (D-005, AC4, chosen at the plan gate).
- blame-history #2: the skip applies below 19 requested surrogates, unlike phase — reject, planned change (AC4 states it, chosen at the plan gate).
- blame-history #3: two-lag pins re-baselined and some p-values moved — reject, planned change (T4).
- blame-history #4: circular cells with different lags no longer share draws — reject, planned change (one matrix per lag, chosen at the plan gate).
- blame-history #5: a circular cell with 19 <= k < `n_surrogates` uses fewer surrogates without a warning — reject, planned change (AC3). The count shows in the `n_surrogates` column and `print()`, as for segment cells.
- blame-history #6: the one-window tests no longer run circular surrogates — reject, planned change (T4 moved them; the window count is still tested).
- blame-history #7: the M018 n = 20 short-series case moved to n = 30 — reject, false as a loss: the n = 20 case is still tested near `test-multiverse.R:903`.
- blame-history #8: the two-matrices test mocks the internal helper — reject, planned change (T4), and AC2 covers per-lag matrices.
- blame-history #9: the autotune message and test moved to the new rule with no coverage loss — reject, no defect reported.
- prior-review #1: the `n_surrogates` grid-column docs contradict NEWS (M017 lesson on counts in prose) — fix now, fixed f8bd33a.
- prior-review #2: the comment above `multiverse_surrogate_text()` names only segment cells — fix now, fixed f8bd33a.
- prior-review #3: `m` and `n - m` in the warning print raw — reject, false: both are integer, so no scientific notation.
- prior-review #4 to #7: checks of cli plural forms, braces, `\code{}`, the print label, and the crossing docs found no problem — reject, no defect reported.

Review pass 2 on 2026-10-11, after review return 1 (T6), at 0623514, main e3d98b9 unmoved.

- Evidence rerun: `test-multiverse-circular.R` 793 expectations, `test-multiverse-pins.R` 11, `test-surrogate.R` 167, all 0 failed. AC1 to AC5 hold on the fixed tree. AC2 holds with the floor, because a shift of at least 1 meets `2 * lag_max <= s` at lag 0. AC6 was re-read after the pass-2 doc fixes (753a620): both Rd files state the shift range with the floor, k, the count, and the skip, and both skip messages state the 19-shift rule.
- Consistency gate (pass 2): `check()` gave 0 errors, 0 warnings, 0 notes. After 753a620, `document()` was current, the full suite gave 3481 expectations with 0 failed, and spelling was clean.

spawned: diff-bug, blame-history, prior-review

- diff-bug #1: the `autotune_wcc()` docs omit the floor of 1 — fix now, fixed 753a620.
- diff-bug #2: the NEWS "same results as before" claim is not exact for a grid whose only lag is 0 — fix now, narrowed to a `lag_max` of at least 1, fixed 753a620.
- diff-bug #3: a lag-0 circular cell draws surrogates but never gives a result — reject, false as a defect: the cell is NA through the estimator's own short-window path, as a phase cell of a 1-sample window is.
- diff-bug #4: a wgranger cell with a 1-sample window shows `n_surrogates` 19 and p `NA` (also for phase) — follow-up, pre-existing, new candidate row "wgranger cells with no p-value report surrogates".
- diff-bug #5: `window_sec = Inf`, or a window past the integer limit, stops the call with "missing value where TRUE/FALSE needed" — follow-up, pre-existing, new candidate row "Non-finite or huge window_sec in the multiverse".
- diff-bug #6: lines past 80 columns in roxygen and NEWS — reject, style (rewrapped while fixing #1 and #2).
- diff-bug #7: no test of the zero-shift message — fix now, test added, fixed 753a620.
- diff-bug #8: no test of a wgranger 1-sample window on the circular path — reject, false as a gap: the AC1 test checks each wgranger shift range from the window, and m = window holds for every window.
- blame-history #1: the floor departs from D-005's `2 * lag_max <= s` — reject, false: shifts of at least 1 meet that bound at lag 0. Recorded in Decisions below.
- blame-history #2: the NEWS range bullet omits the floor — fix now, fixed 753a620.
- blame-history #3: the M018 `n_surrogates` NEWS bullet now also names circular cells — reject, false as a defect: the development NEWS describes the next release as a whole, and the text is true of it.
- blame-history #4 to #6: line length and ragged wraps — reject, style.
- blame-history #7: the zero-shift message branch undoes nothing — reject, no defect reported.
- blame-history #8: the lag-0 test guards a defect that only this branch had — reject, false as a defect: the comment says "Without the floor", which is accurate.
- prior-review #1: the autotune docs keep the unfloored rule (M017 lesson on counts in prose) — fix now, fixed 753a620.
- prior-review #2: the autotune skip message states no floor — fix now, fixed 753a620.
- prior-review #3: the NEWS same-results claim and lag 0 — fix now, same fix as diff-bug #2, fixed 753a620.
- prior-review #4 and #5: line length — reject, style.

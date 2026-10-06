<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section. -->
# M012: One significance rule (p <= .05) and surrogate print edges

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Resolves:** —
- **Surface tier:** user-facing — changes the significance calls that print methods, `synchrony_multiverse()`, and `autotune_wcc()` report
- **Branch/PR:** m012-significance-rule

## Goal

Make every surrogate significance call in bsync come from one internal p <= .05 rule that also handles NA p-values and unreachable thresholds.

## Scope

**In:**
- One internal helper, `is_significant()`, returns `!is.na(p) & p <= 0.05`. The add-one p-value controls the size of a test that rejects at p <= alpha (phipson2010, pp. 4, 6, and D-002).
- These surrogate significance calls go through the helper: the four surrogate print methods (`R/surrogate_analysis.R:619,656,695,712,749`) and the multiverse summary (`R/multiverse.R:398,412`). So do `plot.bsync_multiverse()` (`R/multiverse_plot.R:46`) and the `autotune_wcc()` significance rate (`R/autotune.R:338`).
- `warn_few_surrogates()` (`R/multiverse.R:23`) moves its threshold from 20 to 19, because 1 / (19 + 1) = .05 now passes.
- The four print methods: an NA p-value gives a message that no significance call was made, not an error. With `n_surrogates` of 18 or fewer, they add a note that no result can reach p <= .05.
- Text: labels, roxygen, `man/`, `vignettes/surrogate-testing.Rmd`, `vignettes/choosing-parameters.Rmd`, and the NEWS.md development section state `p <= .05`, the 19-surrogate floor, and the 20-dyad floor for `keep_roles = TRUE` (11 dyads for `keep_roles = FALSE` is unchanged).

**Out:**
- Per-window parametric Granger p-values in `wgranger()` print and plot (`R/wgranger.R:96-100`, `R/wgranger_plot.R:78`) keep `< 0.05`. They are continuous F-test p-values, not add-one surrogate p-values. For them, p = .05 exactly has probability zero. No row.
- A call-time warning in the four `*_surrogate()` wrappers: rejected at plan (work log). The print note gives the message.
- An `alpha` argument: not asked for, no row. The helper leaves one site to change.

## Acceptance criteria

- [x] AC1: A unit test shows `is_significant(c((4 + 1) / (99 + 1), 1 / 20, 0.0501, NA_real_, 0.01))` returns `c(TRUE, TRUE, FALSE, FALSE, TRUE)`.
- [x] AC2: Four text sweeps hold.
      (a) `grep -rnE '<=? *0\.05' R/` lists exactly one line outside `R/wgranger.R` and `R/wgranger_plot.R`, and that line is in the body of `is_significant()`.
      (b) `grep -rnE '< *0?\.05' R/ man/ vignettes/ NEWS.md README.Rmd` lists only lines in `R/wgranger.R`, `R/wgranger_plot.R`, and `man/` pages generated from them.
      (c) `grep -rnE '21 dyads|[Bb]elow 20|n_surrogates *[<>]=? *20' R/ man/ vignettes/ NEWS.md` lists nothing.
      (d) `grep -n 'p <= \.05'` lists at least one line in each of `R/multiverse.R` and `R/multiverse_plot.R`. `grep -n '20 dyads'` lists at least one line in each of `R/surrogate_generation.R`, `vignettes/surrogate-testing.Rmd`, and `NEWS.md`.
- [x] AC3: Tests show that `synchrony_multiverse()` warns with class `bsync_few_surrogates` at `n_surrogates = 18` and not at 19. Tests show that `autotune_wcc()` gives that warning once at 18 and not at 19.
- [x] AC4: Tests show that a p-value of exactly .05 counts as significant, and an NA p-value does not, in three places: `n_significant` of the multiverse summary, the significance marking of `plot.bsync_multiverse()`, and the significance rate of `select_specification()`.
- [x] AC5: Tests cover six NA cases: `print.wcc_surr()`, `print.wdtw_surr()`, and `print.wphase_surr()` with an NA p-value, and `print.wgranger_surr()` with `p_value_xy` NA, with `p_value_yx` NA, and with both NA. In each case, printing gives no error and returns `x` invisibly. Each NA p-value gives a message that no significance call was made. For wgranger, that message names the NA direction, and the other direction's call still prints.
- [x] AC6: Tests show that each print method calls p = .05 significant, wgranger in both directions. Tests also show that each print method prints the new note that no result can reach p <= .05 at `n_surrogates = 18`, and not at 19. The tests match the note by its own text.
- [ ] AC7: `devtools::test()` passes, `devtools::document()` leaves no diff, `devtools::check()` gives 0 errors and 0 warnings, and the NEWS.md development section states the `p <= .05` rule.

## Coverage

- AC1 → T1
- AC2 → T2, T3, T4
- AC3 → T3
- AC4 → T3
- AC5 → T2, T6
- AC6 → T2
- AC7 → T4, T5

## Tasks

- [x] T1: Add `is_significant()` in `R/surrogate_analysis.R` beside `format_p_value()`, with a comment citing phipson2010 and D-002, and its unit test (AC1). Keep `0.05` out of the comment so sweep (a) finds one line.
- [x] T2: Route the five comparisons in the four print methods through `is_significant()`. Add the NA-p message and the few-surrogates note to each. Extend the `mock_surr()` tests in `tests/testthat/test-surrogate.R` per class, asserting with `expect_message` (cli writes to the message stream). The existing `< 1000` note prints at 18 and at 19, so match the new note by its own text.
- [x] T3: Route `R/multiverse.R:398,412`, `R/multiverse_plot.R:46`, and `R/autotune.R:338` through the helper, with the AC4 tests. Move `warn_few_surrogates()` (`R/multiverse.R:24,32`) to 19 and update its message and comment. Update the labels at `R/multiverse.R:566` and `R/multiverse_plot.R:63`. Rewrite the test at `tests/testthat/test-multiverse.R:52` for 18 and 19, and add the `autotune_wcc()` case. Re-record the vdiffr snapshot `_snaps/multiverse/multiverse-spec-curve.svg` and make sure that its diff changes only the legend label.
- [x] T4: Update roxygen (`R/multiverse.R:87`, `R/autotune.R:83-90,289`, `R/multiverse_plot.R:14,18`, `R/surrogate_generation.R:177-187`), then run `devtools::document()`. Update `vignettes/surrogate-testing.Rmd` (lines 45, 117, 123, 153, including the code `sum(p_circular < 0.05)`) and `vignettes/choosing-parameters.Rmd:214`. Rewrite the NEWS.md development entries in place. State the new rule without quoting the old `<` form, because sweep (b) covers NEWS.md.
- [x] T5: Run `spelling::spell_check_package()`, `devtools::test()`, and `devtools::check()`. Run `air format` only on touched files that were air-clean before the edit.
- [x] T6: Review return 1: add `expect_invisible()` for the `p_value_xy`-NA and `p_value_yx`-NA wgranger cases in the NA print test (AC5).

## Work log

- 2026-10-06: created by /milestone-plan.
- 2026-10-06: question set: which work to plan — the significance-edges cluster (the two M011 review candidate rows).
- 2026-10-06: question set: does p = .05 exactly count as significant — yes, use p <= .05.
- 2026-10-06: collision check: absorbs candidate rows "Significance boundary for add-one p-values" and "Surrogate print methods at the edges". D-001's "at least 21 dyads" becomes 20 under the new rule, so D-002 narrows it. No open GitHub issues or PRs.
- 2026-10-06: plan gate chose a print-time note over a call-time warning in the four wrappers because a call-time warning fires once per dyad in a pseudo-dyad loop; falsified by a user report of an unreachable p-value read from a wrapper without printing it.
- 2026-10-06: lessons surfaced: cli output asserts with expect_message (M008). `R/autotune.R` and `tests/testthat/test-autotune.R` are not air-clean (M010). New words need `inst/WORDLIST` (M011). Vignettes need `plan(sequential)` under `load_all()` (M011).
- 2026-10-06: criteria audit (full mode, fresh Opus reader) returned 11 findings, all fixed. AC1 adds a vector probe. Old AC2 and AC3 merge into one sweep AC with a one-line helper rule, a wider threshold pattern, and positive checks for the new text. New AC4 tests the multiverse, plot, and `select_specification()` counts. AC3 adds `autotune_wcc()`. AC5 adds the both-NA wgranger case and direction naming. AC6 tests both wgranger directions and matches the new note by text. T3 re-records the vdiffr snapshot. T4 names the vignette code lines and NEWS wording. D-002 is written in the plan commit. The audit found no Invariant conflict and no float issue at p = .05.
- 2026-10-06: T1 done: `is_significant()` and `significance_reachable()` added in `R/surrogate_analysis.R`, with a unit test. Planting the old `<` rule makes the AC1 assertion fail. `R/surrogate_analysis.R` and `tests/testthat/test-surrogate.R` were not air-clean before the edit, so air is not run on them.
- 2026-10-06: T2 done: the four print methods call shared helpers `print_significance_call()` and `print_surrogate_count_notes()` (implementation choice: one place for the NA message and both notes). The NA, p = .05, and 18/19-note tests fail against the old print code (1 error, 5 and 4 failures) and pass on the branch.
- 2026-10-06: T3 done: the multiverse summary moved into `multiverse_robustness()` and the plot's data step into `multiverse_plot_data()` (implementation choice: internal functions the AC4 tests call directly). `warn_few_surrogates()` uses `significance_reachable()`. With the old `<` rule planted, the 18/19, AC4, and autotune tests fail. The re-recorded snapshot changes only the legend text, and its 14 points keep their colors (7 and 7).
- 2026-10-06: T4 done: roxygen, five `man/` pages, both vignettes, and the NEWS.md development section now state p <= .05, 19 surrogates, and 20 dyads. The four AC2 sweeps give the planned output: (a) and (b) list only the helper body line and `R/wgranger.R:96,97,100`, (c) lists nothing, and (d) finds each required line.
- 2026-10-06: T5 in progress: spelling clean after the vignette's `$p \le .05$` became code, and `devtools::test()` passes with no failures or skips. `devtools::check()` is running. Checkpoint before the claim-audit re-read and the check result.
- 2026-10-06: claim audit: 56 claims read, 4 corrected — tests/testthat/test-surrogate.R, R/autotune.R (+ man/autotune_wcc.Rd), NEWS.md, cairn/references/phipson2010.md. The same pass corrected the stale "below 20" in cairn/DESIGN.md:247. Re-read of the corrections pending.
- 2026-10-06: claim audit re-read: all 5 corrections OK. NEWS now says "computable cell", per the reader's optional point.
- 2026-10-06: T5 done: `devtools::check()` gave 0 errors, 0 warnings, and 0 notes, on the tree before the claim-audit text fixes. After those fixes, `devtools::document()` leaves no diff and `tools::checkRd()` is clean on `man/autotune_wcc.Rd`. Air ran only on `R/surrogate_generation.R`, the one touched file that was air-clean before. Status set to review.
- 2026-10-06: review return 1: AC5 fails as written. The `p_value_xy`-NA and `p_value_yx`-NA wgranger cases in "print methods make no significance call on an NA p-value" (`tests/testthat/test-surrogate.R`) never assert that `print()` returns `x` invisibly.
- 2026-10-06: T6 added for review return 1 (minor amendment: new task, Coverage AC5 → T2, T6). T6 done: both single-NA wgranger cases now assert `expect_no_error` and `expect_invisible`. The test passes with 27 expectations. Status set to review.

## Decisions

## Review

Pass 1, 2026-10-06, branch head 95b00b0, synced with origin/main 9eff39c.

- AC1: `test-surrogate.R` "is_significant() counts p <= .05 and never NA" passes (2 expectations). Its first expectation is the AC1 vector, word for word. With the old `<` rule planted, it fails (T1 work log).
- AC2: (a) lists `R/surrogate_analysis.R:592`, the body line of `is_significant()`, and `R/wgranger.R:96,97,100`. (b) lists only `R/wgranger.R:96,97,100`. (c) lists nothing (grep exit 1). (d) `p <= .05` count: `R/multiverse.R` 5, `R/multiverse_plot.R` 3. `20 dyads` count: `R/surrogate_generation.R` 1, `vignettes/surrogate-testing.Rmd` 2, `NEWS.md` 1.
- AC3: `test-multiverse.R` "fewer than 19 surrogates warns ..." passes: a `bsync_few_surrogates` warning at 18 and `expect_no_condition(class =)` at 19. `test-autotune.R` "autotune_wcc gives the few-surrogates warning once per call" passes: 1 warning at 18 and 0 at 19, with `NOT_CRAN=true`.
- AC4: `test-multiverse.R` "the multiverse counts p = .05 as significant and NA p as not" passes. `n_significant` is 1 of 3 valid cells (p = .05 counted, NA p not counted), and `multiverse_plot_data()` marks `c(FALSE, TRUE)` for the NA and .05 cells. `test-autotune.R` "select_specification counts p = .05 ..." passes: the rate is 0.5 for p = .05, NA, .5, so the NA dyad is not counted as significant.
- AC5: FAILS as written. "print methods make no significance call on an NA p-value" passes, but it asserts `expect_invisible()` only for wcc, wdtw, wphase, and wgranger with both NA. The `p_value_xy`-NA and `p_value_yx`-NA wgranger cases never check that `print()` returns `x` invisibly. Returned to implement (review return 1).
- AC6: "print methods call p = .05 significant" passes (wcc, wdtw, wphase, and wgranger in both directions). "print methods note when no result can reach p <= .05" passes for all four classes at 18 (present) and 19 (absent), matched by the note's own text "so no result can reach p <= .05".

Pass 2, 2026-10-06, branch head 018c0ec, still synced with origin/main 9eff39c.

- AC5: "print methods make no significance call on an NA p-value" passes with 27 expectations. It covers wcc, wdtw, and wphase with NA p, and wgranger with `p_value_xy` NA, with `p_value_yx` NA, and with both NA. Each case asserts no error, the "No significance call" message, and `expect_invisible()`. The wgranger cases also assert the named direction and the other direction's call.
- Consistency gate: `cairn_validate.py` exit 0 (all checks passed). `devtools::document()` leaves no diff. README.Rmd is untouched by the branch. `pkgdown::check_pkgdown()` reports no problems. The branch adds no top-level files. NEWS.md has the development-section entries. `NOT_CRAN=true devtools::test()`: 17 files, 243 tests, 0 failed, 0 errors, 0 skipped. `devtools::check()` at 018c0ec: 0 errors, 0 warnings, 0 notes.
- spawned: diff-bug, blame-history, prior-review
- diff-bug #1: `plot.bsync_multiverse()` now draws an NA-p cell as not significant (before, NA color), and NEWS does not say so — fix now (NEWS clause).
- diff-bug #2: `print_significance_call()` passes `yes`/`no` to cli as format strings, so a brace in a future message errors (probe confirmed "object 'b' not found") — fix now (`"{yes}"`/`"{no}"`, plus a test).
- diff-bug #3: an NA `n_surrogates` prints a misleading "With NA surrogate" note before the old `< 1000` error. Only hand-built objects reach it — fix now (skip the note on NA).
- diff-bug #4: a vector `p` errors in `print_significance_call()` — reject, false as a defect: every surrogate object holds a scalar p, and the old code behaved the same.
- diff-bug #5: the autotune 19-surrogate case asserts absence by message text, and the AC4 plot test calls `multiverse_plot_data()` on a bare list — reject, false as a defect: the 19 case collects every warning the call gives, and `plot.bsync_multiverse()` calls `multiverse_plot_data()` (`R/multiverse_plot.R`), which the snapshot test runs.
- diff-bug #6: `R/surrogate_generation.R:183` is an 87-character roxygen line — reject, style.
- diff-bug #7: for `estimator = "wgranger"`, the multiverse summary and plot read only the x -> y `p` and ignore `p_yx`. The code before the branch did the same — follow-up, new candidate row "Multiverse summary and plot ignore the y -> x Granger direction".
- blame-history #1: the re-recorded snapshot changes legend x-coordinates and rewrites tick labels — reject, false: the 32 SVG text labels match except the legend title, and the tick labels (2.5, 10.0, 12.5) are present. Coordinates move because the wider legend title narrows the panel by 0.65 px.
- blame-history #2: the plot NA-p color change, which NEWS does not mention — fix now (same fix as diff-bug #1).
- blame-history #3: the new note prints before the `< 1000` note, changing output order — reject, false as a defect: no test or doc depends on the order, and the suite is green.
- blame-history #4: `format_p_value()` rounding survives, and the call uses the raw p — reject, false: the reviewer reports no defect.
- blame-history #5: `wgranger()` per-window p-values keep `< 0.05` — reject, planned change (Scope Out, D-002).
- blame-history #6: `sign_consistent`, the `n_cells`/`n_valid` meanings, the warning class, and the autotune muffle are unchanged — reject, false: no defect reported.
- blame-history #7: the 19 threshold has no floating-point edge — reject, false: no defect reported.
- blame-history #8: the long roxygen line, D-001's old text, and unticked AC5/AC7 boxes — reject: the line is style, D-001 is append-only history narrowed by D-002, and the boxes are expected at the gate.
- prior-review #1: the wrappers still give no call-time note for too few surrogates — reject, planned change (plan-gate choice, work log).
- prior-review #2: the per-window Granger summary keeps `< 0.05` — reject, planned change (Scope Out, D-002).

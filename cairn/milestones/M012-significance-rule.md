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

- [ ] AC1: A unit test shows `is_significant(c((4 + 1) / (99 + 1), 1 / 20, 0.0501, NA_real_, 0.01))` returns `c(TRUE, TRUE, FALSE, FALSE, TRUE)`.
- [ ] AC2: Four text sweeps hold.
      (a) `grep -rnE '<=? *0\.05' R/` lists exactly one line outside `R/wgranger.R` and `R/wgranger_plot.R`, and that line is in the body of `is_significant()`.
      (b) `grep -rnE '< *0?\.05' R/ man/ vignettes/ NEWS.md README.Rmd` lists only lines in `R/wgranger.R`, `R/wgranger_plot.R`, and `man/` pages generated from them.
      (c) `grep -rnE '21 dyads|[Bb]elow 20|n_surrogates *[<>]=? *20' R/ man/ vignettes/ NEWS.md` lists nothing.
      (d) `grep -n 'p <= \.05'` lists at least one line in each of `R/multiverse.R` and `R/multiverse_plot.R`. `grep -n '20 dyads'` lists at least one line in each of `R/surrogate_generation.R`, `vignettes/surrogate-testing.Rmd`, and `NEWS.md`.
- [ ] AC3: Tests show that `synchrony_multiverse()` warns with class `bsync_few_surrogates` at `n_surrogates = 18` and not at 19. Tests show that `autotune_wcc()` gives that warning once at 18 and not at 19.
- [ ] AC4: Tests show that a p-value of exactly .05 counts as significant, and an NA p-value does not, in three places: `n_significant` of the multiverse summary, the significance marking of `plot.bsync_multiverse()`, and the significance rate of `select_specification()`.
- [ ] AC5: Tests cover six NA cases: `print.wcc_surr()`, `print.wdtw_surr()`, and `print.wphase_surr()` with an NA p-value, and `print.wgranger_surr()` with `p_value_xy` NA, with `p_value_yx` NA, and with both NA. In each case, printing gives no error and returns `x` invisibly. Each NA p-value gives a message that no significance call was made. For wgranger, that message names the NA direction, and the other direction's call still prints.
- [ ] AC6: Tests show that each print method calls p = .05 significant, wgranger in both directions. Tests also show that each print method prints the new note that no result can reach p <= .05 at `n_surrogates = 18`, and not at 19. The tests match the note by its own text.
- [ ] AC7: `devtools::test()` passes, `devtools::document()` leaves no diff, `devtools::check()` gives 0 errors and 0 warnings, and the NEWS.md development section states the `p <= .05` rule.

## Coverage

- AC1 → T1
- AC2 → T2, T3, T4
- AC3 → T3
- AC4 → T3
- AC5 → T2
- AC6 → T2
- AC7 → T4, T5

## Tasks

- [x] T1: Add `is_significant()` in `R/surrogate_analysis.R` beside `format_p_value()`, with a comment citing phipson2010 and D-002, and its unit test (AC1). Keep `0.05` out of the comment so sweep (a) finds one line.
- [x] T2: Route the five comparisons in the four print methods through `is_significant()`. Add the NA-p message and the few-surrogates note to each. Extend the `mock_surr()` tests in `tests/testthat/test-surrogate.R` per class, asserting with `expect_message` (cli writes to the message stream). The existing `< 1000` note prints at 18 and at 19, so match the new note by its own text.
- [x] T3: Route `R/multiverse.R:398,412`, `R/multiverse_plot.R:46`, and `R/autotune.R:338` through the helper, with the AC4 tests. Move `warn_few_surrogates()` (`R/multiverse.R:24,32`) to 19 and update its message and comment. Update the labels at `R/multiverse.R:566` and `R/multiverse_plot.R:63`. Rewrite the test at `tests/testthat/test-multiverse.R:52` for 18 and 19, and add the `autotune_wcc()` case. Re-record the vdiffr snapshot `_snaps/multiverse/multiverse-spec-curve.svg` and make sure that its diff changes only the legend label.
- [x] T4: Update roxygen (`R/multiverse.R:87`, `R/autotune.R:83-90,289`, `R/multiverse_plot.R:14,18`, `R/surrogate_generation.R:177-187`), then run `devtools::document()`. Update `vignettes/surrogate-testing.Rmd` (lines 45, 117, 123, 153, including the code `sum(p_circular < 0.05)`) and `vignettes/choosing-parameters.Rmd:214`. Rewrite the NEWS.md development entries in place. State the new rule without quoting the old `<` form, because sweep (b) covers NEWS.md.
- [x] T5: Run `spelling::spell_check_package()`, `devtools::test()`, and `devtools::check()`. Run `air format` only on touched files that were air-clean before the edit.

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

## Decisions

## Review

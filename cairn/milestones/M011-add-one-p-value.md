<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. The one size check that can fail is
     cairn_validate's <150 over the plan-owned body. -->
# M011: Add-one surrogate p-values

- **Status:** review   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate; M<xx>, M<yy> or — -->
- **Driving RR:** —   <!-- owner: plan · create/amend-via-gate -->
- **Principles touched:** —   <!-- owner: plan · create/amend-via-gate; bsync's numbered invariants live in CLAUDE.md (2 and 3 bind here), no DESIGN IP/GP numbering yet -->
- **Resolves:** —   <!-- owner: plan · create/amend-via-gate -->
- **Surface tier:** user-facing — changes the p-value that every exported surrogate wrapper returns   <!-- owner: plan · create/amend-via-gate -->
- **Branch/PR:** m011-add-one-p-value   <!-- owner: implement (branch) / review (PR URL) · create -->

## Goal
<!-- owner: plan · create; a wrong goal returns to plan, never edited in place -->

Make every surrogate test report the add-one p-value (b + 1) / (n + 1)
(phipson2010 (p. 6)), so that a reported p-value is never 0 and the test
keeps its nominal size.

## Scope
<!-- owner: plan · create/amend-via-gate -->

**In:** The p-value sites in `R/surrogate_analysis.R` (lines 116, 274, 385,
386, 514 at plan time). The `p_value == 0` display branches of the four
print methods (lines 540, 595, 636, 641, 700). A cli warning in
`synchrony_multiverse()` when `n_surrogates < 20`, because the multiverse
and `autotune_wcc()` call a cell significant at p < .05
(`R/multiverse.R:374`, `R/autotune.R:308`). The phipson2010 source note.
DESIGN.md §6 step 3. Roxygen, vignettes, README, and NEWS text that states
how the p-value is computed or its smallest value. D-001 records the
decision (plan commit).

**Out:** The fixed .05 threshold in the multiverse and `autotune_wcc()`
stays as it is. Within-dyad surrogate generators are not changed. IAAFT
and segment shuffling stay on their candidate row.

## Acceptance criteria
<!-- owner: plan · create/amend-via-gate; review reads, never reinterprets. -->

- [x] AC1: `wcc_surrogate()`, `wdtw_surrogate()`, `wgranger_surrogate()`
      (`p_value_xy` and `p_value_yx`), and `wphase_surrogate()` return
      p = (b + 1) / (n + 1), phipson2010 (p. 6). Here n is the number of
      surrogate columns, and b is the number of surrogate statistics at least
      as extreme as the observed statistic, in the tail each wrapper uses
      today: upper for wcc, wgranger, and wphase, lower for wdtw. For each
      wrapper, a test counts b with explicit code from the returned observed
      and surrogate statistics and asserts the p-value to 1e-12. The claim
      covers inputs where every observed and surrogate statistic is finite.
      NA handling is unchanged. The probes include one case with b = 0 and
      one with 0 < b < n for each wrapper, and for both wgranger directions.
- [x] AC2: A simulation test shows that the add-one p-value keeps its size
      under an exchangeable null, phipson2010 (pp. 4, 6). Each replicate
      builds N = 10 independent, uncoupled dyads of white noise and tests
      dyad 1 against all 9 of its `generate_surrogate_pseudo()` partners
      (`keep_roles = TRUE`) with `wcc_surrogate()`. The test runs at least
      200 replicates under a fixed seed. Under this null the p-value takes
      the values k / 10 with equal chance. The test asserts that no p-value
      is below 0.1, and that the shares of p-values at or below 0.15 and at
      or below 0.5 are within 3 binomial standard errors of 0.1 and 0.5.
      (Under b / n the share at or below 0.15 is 0.2.)
- [x] AC3: Each of the four print methods shows the returned p-value
      rounded to 4 digits, or `< 0.0001` when that rounding gives 0. The
      `< 1/n` form is gone. Tests on constructed result objects assert the
      printed p-value (captured from the message stream) for a mid-range p
      and for a p below 5e-5, for each print method and for both wgranger
      directions.
- [x] AC4: `synchrony_multiverse()` raises a cli warning when
      `n_surrogates < 20`. The warning states that the smallest possible
      p-value is `1 / (n_surrogates + 1)`, so no cell can reach p < .05.
      `autotune_wcc()` gives this warning once per call. Tests assert the
      warning with a message matcher at `n_surrogates = 19` and no such
      warning at `n_surrogates = 20`, for both functions, and assert that
      `autotune_wcc()` with 3 dyads and 19 surrogates gives it once.
- [x] AC5: `cairn/references/phipson2010.md` exists with its `INDEX.md`
      line and anchors the formula (p. 6) and the size of b / n (p. 4).
      DESIGN.md §6 states the add-one formula with that citation in step 3
      and in its Invariant 2 paragraph. The roxygen of the four wrappers
      states the formula and why it never gives 0. The surrogate-testing
      vignette and the `generate_surrogate_pseudo()` roxygen state that with
      N dyads and `keep_roles = TRUE` the smallest single-dyad pseudo-dyad
      p-value is 1 / N, so p < .05 needs at least 21 dyads. The wcc and wdtw workflow vignettes no
      longer report a p-value as `< .001` or `< 0.001`. NEWS.md has an
      entry. README.md is re-knit from README.Rmd.
- [x] AC6: `devtools::document()` leaves no diff, `devtools::test()`
      passes, and `devtools::check()` reports 0 errors and 0 warnings.

## Coverage
<!-- owner: plan · create/amend-via-gate -->

- AC1 → T1
- AC2 → T1
- AC3 → T2
- AC4 → T2
- AC5 → T3, T4
- AC6 → T4

## Tasks
<!-- owner: plan (create) / implement (check-off, minor edits). -->

- [x] T1: Tests first in `tests/testthat/test-surrogate.R` (AC1) and
      `test-surrogate-pseudo.R` (AC2). Change the five p-value lines in
      `R/surrogate_analysis.R` to `(b + 1) / (n_surrogates + 1)`. Update the
      seeded pin (`test-surrogate.R:462`) to the add-one values from the same
      seed, with `(b + 1) / (99 + 1)` in a comment. Fix any other test whose
      expected p-value the change moves.
- [x] T2: Tests first. Replace the `p_value == 0` branches of the print
      methods with the AC3 display. Add the `n_surrogates < 20` warning to
      `synchrony_multiverse()`, and make `autotune_wcc()` give it once (for
      example, warn up front and muffle the per-dyad copies). Update tests in `test-multiverse.R` and
      `test-autotune.R` that use fewer than 20 surrogates: raise the count or
      expect the warning, and keep any assertion that needs a significant
      cell valid.
- [x] T3: Author `cairn/references/phipson2010.md` from the source-note
      template, with its `INDEX.md` line (PDF on the shelf as
      `sources/phipson2010.pdf`). Update DESIGN.md §6 (step 3 and the
      Invariant 2 paragraph).
- [x] T4: As a finding aid, run `grep -rniE 'p-value|p_value|p value|tail'`
      over `R/`, `vignettes/`, and `README.Rmd`, and log the hit count and
      the hits changed. Update the AC5 sites: wrapper roxygen, the
      `generate_surrogate_pseudo()` "How many surrogates" paragraph, the
      surrogate-testing vignette (pseudo-dyad demo and its text), and the
      wcc and wdtw workflow prose that says "< .001". Fix stale code
      comments found on the way (`R/multiverse.R:464`). Re-knit README. Add
      NEWS. Run `devtools::document()`, `devtools::test()`,
      `devtools::check()`, and `air format`.

## Work log
<!-- owner: any skill · append-only; one line per entry; absolute dates. -->

- 2026-10-05: created by /milestone-plan.
- 2026-10-05: question set: switch surrogate p-values from b/n to (b+1)/(n+1) as a default change — "Yes, add to this plan (Recommended)". Recorded as D-001.
- 2026-10-05: lineage: absorbs the candidate row "Add-one surrogate p-value" (M009 review); row removed in the plan commit.
- 2026-10-05: source: Phipson & Smyth 2010, arXiv:1603.05766 (open access, 199.6 KB), fetched 2026-10-05 and shelved as `cairn/references/sources/phipson2010.pdf`. Formula read on PDF p. 6 (§4), size of b/m on p. 4 (§2).
- 2026-10-05: plan gate chose a warning in `synchrony_multiverse()` for `n_surrogates < 20` over no warning because without it `autotune_wcc()` fails its gate with only the generic "no specification passed" message; falsified by users who run small surrogate counts on purpose and find the warning noise.
- 2026-10-05: plan gate chose the same add-one form for the exhaustive pseudo-dyad null over a separate rank form because with every partner used, (b + 1) / (n + 1) is the observed partner's rank among the N series, so one formula covers both; falsified by a source that defines a different exact p-value for exhaustive enumeration.
- 2026-10-05: criteria audit (full mode, fresh Opus reader) returned 6 findings on M011, all fixed. Print display gives `< 0.0001` where rounding to 4 digits would print 0 (AC3). The simulation tests one dyad per replicate, so the p-values are independent, with at least 200 replicates and a 0.15 cutoff that separates b/n from the add-one form (AC2). AC1 is scoped to finite statistics and probes b = 0 in both wgranger directions. `autotune_wcc()` warns once (AC4). The grep-sweep domain and the work-log recording act moved from the criterion to T4, and AC5 now names its doc sites, including DESIGN §6's Invariant 2 paragraph. The seeded pin criterion (an instrument) folded into T1.
- 2026-10-05: plan gate left CLAUDE.md Invariant 2 unedited because its wording ("the p-value is its tail") still holds under the add-one form; falsified by a reviewer reading it as b / n.
- 2026-10-05: implement started. Branch m011-add-one-p-value cut from main at f06719a, in sync with origin.
- 2026-10-05: T1 done. The five p-value lines now use (b + 1) / (n + 1). AC1 tests count b with a loop for all five statistics, with 19 circular surrogates. A bidirectional VAR(1) dyad gives b = 0 for all five, and an independent dyad gives 0 < b < 19 for all five. The AC2 test uses 400 replicates (about 10 s, skip_on_cran). Run against b / n, it failed the p >= 0.1 check and the 0.15 share (0.205). The seeded pin is now (91 + 1) / 100 and (83 + 1) / 100. The multiverse snapshot test moved from 9 to 20 surrogates because at 9 no cell can reach p < .05. The new snapshot shows 6 of 12 cells significant, median ES 1.41.
- 2026-10-05: T2 done. The four print methods use a new internal `format_p_value()`, which rounds to 4 digits in fixed notation and shows `< 0.0001` when that rounding gives 0. The test also pins p = 1e-4 as "0.0001" because `as.character()` printed "1e-04". The new `warn_few_surrogates()` in `R/multiverse.R` gives a cli warning of class `bsync_few_surrogates`. `autotune_wcc()` gives it once up front and muffles that class in the per-dyad calls. With the class removed as a planted defect, the once-per-call test saw 4 warnings. Multiverse tests at 9 or 19 surrogates moved to 20, the odd-length phase test moved from 10 to 20, and the 5-surrogate autotune test now expects the warning. Full suite clean (NOT_CRAN).
- 2026-10-05: T1 and T2 share one checkpoint commit because T2 edits started before the T1 commit. The T4 wrapper roxygen edits in `R/surrogate_analysis.R` ride along in that commit. They add `@md` to the four wrapper blocks so the formula and the journal title render.
- 2026-10-05: T3 done. `cairn/references/phipson2010.md` was written from PDF pp. 1 and 3 to 7, read directly. It anchors b / m (p. 3), its size (p. 4), the add-one formula (p. 6), and the exhaustive form (p. 7). `INDEX.md` has its line. DESIGN.md §6 step 3 states the formula, tails, and citation, and the Invariant 2 paragraph names the add-one form.
- 2026-10-05: T4 in progress, checkpoint. The finding-aid grep found 88 hits before the edits and 103 after. Changed: four wrapper roxygen blocks, the `generate_surrogate_pseudo()` paragraph, the `R/multiverse.R` WDTW comment, and the surrogate-testing, wcc, wdtw, and wgranger vignettes. The wcc and pseudo-dyad prose was written against a sequential run of the vignette code: the wcc example gives b = 0, p = 1/1001, printed 0.001, and the circular null flags 7 of 10 dyads. README re-knit shows 0.0693. NEWS entry added, and Phipson and Smyth added to `inst/WORDLIST`. `air format` was not kept because it restyled about 26 untouched files. A candidate row records that. `devtools::check()` is running.
- 2026-10-05: T4 done. `devtools::check()` gave Status OK, with 0 errors, 0 warnings, and 0 notes. A first run had one test-output NOTE because the spelling test ran before the WORDLIST edit. `devtools::document()` leaves no diff.
- 2026-10-05: delegation: a fresh Opus reader (general-purpose agent) ran the claim audit. It read the added lines, ran the vignette code sequentially, and ran the changed test files.
- claim audit: 86 claims read, 4 corrected — vignettes/wgranger-workflow.Rmd, vignettes/wdtw-workflow.Rmd, vignettes/surrogate-testing.Rmd, R/surrogate_generation.R
- 2026-10-05: claim audit fixes. The F-statistics in the wgranger prose are now 159.48 and 133.12. The wdtw hard-coded output and prose are now 27.043 and 57.0134. Both were also stale on main. The size sentence in the surrogate-testing vignette now names the exchangeability condition. The pseudo-dyad roxygen rank sentence now names `keep_roles = TRUE`. The same reader re-read all four once, and all four passed.
- 2026-10-05: the plan gate approved the switch as a default change on a 0.x package, and no argument exists to deprecate. That approval is taken as the pre-1.0 waiver of a deprecation cycle. NEWS states the change.
- 2026-10-05: all tasks done, status review.
- 2026-10-05: review fix-now commit. A surrogate matrix with no columns now aborts in all four wrappers. `format_p_value()` ignores `getOption("digits")`. A missing `n_surrogates` gives a cli abort in both functions. The few-surrogates warning now comes after every argument check. Other changes cover the `n_surrogates` docs, DESIGN §9, the 1 / N floor scope, and the `fast_method` caveat. The muffle and line-wrap tests were hardened. Each new test failed before its fix, and the full suite passes.

## Decisions
<!-- owner: implement / review · append-only; milestone-local. -->

## Review
<!-- owner: review · exclusive -->

Evidence gathered 2026-10-05 at 75f348f (main f06719a, unmoved), NOT_CRAN=true.

- AC1 evidence: `test-surrogate.R` tests "all four wrappers return (b + 1) / (n + 1) when b = 0" (11 expectations) and "... when 0 < b < n" (10 expectations) pass, 0 failed, 0 skipped. Each counts b with a plain loop for wcc, wdtw (lower tail), wgranger xy and yx, and wphase, and checks the p-value to 1e-12. Every observed and surrogate statistic in the probes is finite. The seeded pin test (2 expectations) also passes.
- AC2 evidence: `test-surrogate-pseudo.R` test "the add-one pseudo-dyad p-value keeps its size under the null" passes (3 expectations, 0 failed, not skipped under NOT_CRAN). It runs 400 seeded replicates of 10 white-noise dyads, with dyad 1 tested against all 9 `keep_roles = TRUE` partners by `wcc_surrogate()`. It asserts every p >= 0.1, and that the shares at or below 0.15 and 0.5 sit within 3 binomial SEs of 0.1 and 0.5. Run against b / n during T1, it failed the first two assertions (share 0.205).
- AC3 evidence: `test-surrogate.R` test "print methods show p rounded to 4 digits, or < 0.0001" passes (11 expectations). It builds result objects for wcc, wdtw, and wphase, and reads the printed value from the message stream. 0.123456 prints "0.1235", 1e-5 prints "< 0.0001", and 1e-4 prints "0.0001". For wgranger, both direction orders are checked. `grep -n '1 / x\$n_surrogates' R/surrogate_analysis.R` finds no match, so the `< 1/n` form is gone.
- AC4 evidence: `test-multiverse.R` test "fewer than 20 surrogates warns that no cell can reach p < .05" passes (2 expectations). At 19 it expects a warning of class `bsync_few_surrogates` whose message matches "smallest possible p-value is `1 / (n_surrogates + 1)` = 0.05". At 20 it expects no such warning. `test-autotune.R` test "autotune_wcc gives the few-surrogates warning once per call" passes (2 expectations). With 3 dyads it counts exactly 1 matching warning at 19 and 0 at 20. With the class removed as a planted defect (T2), the count was 4.
- AC5 evidence: `cairn/references/phipson2010.md` exists, and `INDEX.md` has its line. The note anchors p. 6 (formula) and p. 4 (size of b / m). DESIGN.md §6 step 3 (line 187) states the formula with the phipson2010 p. 6 citation, and the Invariant 2 paragraph (line 195) names the add-one form with that citation. All 4 wrapper roxygen blocks contain "(b + 1) / (n + 1) (Phipson & Smyth, 2010)" and "never 0". The surrogate-testing vignette (lines 45, 123) and the `generate_surrogate_pseudo()` roxygen (line 180) state the 1 / N floor and 21 dyads. A grep for `< .001` or `< 0.001` in the wcc and wdtw workflow vignettes finds 0 matches. NEWS.md line 3 has the entry. README.md was re-knit in 2d9e73d and shows 0.0693.
- AC6 evidence: `devtools::document()` at 75f348f leaves `git status` clean. `devtools::check()` at 75f348f gives Status OK, with 0 errors, 0 warnings, and 0 notes. That check runs the full test suite, so `devtools::test()` passes too.
- Consistency gate: `cairn_validate.py` passes all checks. `pkgdown::check_pkgdown()` reports no problems. README.md is knit from README.Rmd. The NEWS.md entry has no milestone numbers, and the doc-hygiene test passes. No new top-level files were added, so no `.Rbuildignore` change is needed. No DESIGN principle changed, so `cairn_impact` is skipped.
- spawned: diff-bug, blame-history, prior-review
- diff-bug #1: a zero-column surrogate matrix gave p = 1 with no message — fix now, fixed 0b040f8
- diff-bug #2: print methods crash on an NA p-value at `if (x$p_value < 0.05)` — follow-up (already on main before M011), candidate row "Surrogate print methods at the edges"
- diff-bug #3: `n_surrogates = NA` gave a base-R error in `autotune_wcc()` (and in `synchrony_multiverse()` on main) — fix now, fixed 0b040f8
- diff-bug #4: `format_p_value()` followed `getOption("digits")`, so `digits = 3` printed "0.123" — fix now, fixed 0b040f8
- diff-bug #5: the 1 / N floor holds only when all partners are used, and `keep_roles = FALSE` was not stated — fix now, fixed 0b040f8
- diff-bug #6: the `n_surrogates` docs and DESIGN §9 did not mention the new warning — fix now, fixed 0b040f8
- diff-bug #7: the once-per-call test did not show that the muffle leaves other warnings alone — fix now, fixed 0b040f8
- diff-bug #8: the warning-text regex depended on cli line wrapping — fix now, fixed 0b040f8
- diff-bug #9: the multiverse warned before its `lag_sec` check, so a call that aborts warned first — fix now, fixed 0b040f8. In `autotune_wcc()` the warning already follows its own argument checks.
- diff-bug #10: with 10 dyads the vignette contrast holds by construction — reject (planned change: Scope names the pseudo-dyad demo and its text, and the prose says so)
- diff-bug #11: the add-one size claim does not hold for `wdtw_surrogate(fast_method = TRUE)` — fix now, fixed 0b040f8
- diff-bug #12: the wphase calibration test counts `<= .05` while the package uses `< .05` — follow-up, candidate row "Significance boundary for add-one p-values"
- blame-history #1: at n = 99, b = 4 gives p = .05, which is not significant under the strict `< .05` rule, where b / n called it significant — follow-up (Scope keeps the .05 rule unchanged), candidate row "Significance boundary for add-one p-values"
- blame-history #2: the vignette demo lost its empirical contrast — reject (planned change, same as diff-bug #10)
- blame-history #3: the "fewer than 21 dyads" advice did not cover `keep_roles = FALSE` — fix now, fixed 0b040f8 (with diff-bug #5)
- blame-history #4: the wrappers give no note below 20 surrogates, and DESIGN §9 and the param docs were not updated — docs part fix now, fixed 0b040f8. The wrapper note is a follow-up, candidate row "Surrogate print methods at the edges"
- blame-history #5: tests raised to 20 surrogates keep their purpose, and the snapshot still shows 6 of 12 — reject (false: the lens itself found no defect)
- blame-history #6: hand-edited vignette numbers cannot be reproduced because the runs are unseeded — reject (false: all three vignettes call `set.seed(2026)`, and the claim auditor reproduced the numbers)
- blame-history #7: the print branches were replaced and nothing was undone — reject (false: no defect)
- prior-review #1: non-ASCII characters on edited roxygen lines — reject (style, already on main; `Encoding: UTF-8` and check gives 0 notes)
- prior-review #2: `@md` was added to only some blocks — reject (false: this follows the M009 lesson, and the candidate row was updated)
- prior-review #3: the `n_surrogates` param docs did not mention the warning — fix now, fixed 0b040f8 (same as diff-bug #6)
- Return floor: no finding shows an acceptance criterion failing. Diff-bug #1 is an input-validation edge case and was fixed on the branch, so status stays review. Verify after the fixes: `devtools::document()`, then the full `devtools::test()` passes (NOT_CRAN).

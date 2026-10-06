<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. The one size check that can fail is
     cairn_validate's <150 over the plan-owned body. -->
# M011: Add-one surrogate p-values

- **Status:** in-progress   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
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

- [ ] AC1: `wcc_surrogate()`, `wdtw_surrogate()`, `wgranger_surrogate()`
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
- [ ] AC2: A simulation test shows that the add-one p-value keeps its size
      under an exchangeable null, phipson2010 (pp. 4, 6). Each replicate
      builds N = 10 independent, uncoupled dyads of white noise and tests
      dyad 1 against all 9 of its `generate_surrogate_pseudo()` partners
      (`keep_roles = TRUE`) with `wcc_surrogate()`. The test runs at least
      200 replicates under a fixed seed. Under this null the p-value takes
      the values k / 10 with equal chance. The test asserts that no p-value
      is below 0.1, and that the shares of p-values at or below 0.15 and at
      or below 0.5 are within 3 binomial standard errors of 0.1 and 0.5.
      (Under b / n the share at or below 0.15 is 0.2.)
- [ ] AC3: Each of the four print methods shows the returned p-value
      rounded to 4 digits, or `< 0.0001` when that rounding gives 0. The
      `< 1/n` form is gone. Tests on constructed result objects assert the
      printed p-value (captured from the message stream) for a mid-range p
      and for a p below 5e-5, for each print method and for both wgranger
      directions.
- [ ] AC4: `synchrony_multiverse()` raises a cli warning when
      `n_surrogates < 20`. The warning states that the smallest possible
      p-value is `1 / (n_surrogates + 1)`, so no cell can reach p < .05.
      `autotune_wcc()` gives this warning once per call. Tests assert the
      warning with a message matcher at `n_surrogates = 19` and no such
      warning at `n_surrogates = 20`, for both functions, and assert that
      `autotune_wcc()` with 3 dyads and 19 surrogates gives it once.
- [ ] AC5: `cairn/references/phipson2010.md` exists with its `INDEX.md`
      line and anchors the formula (p. 6) and the size of b / n (p. 4).
      DESIGN.md §6 states the add-one formula with that citation in step 3
      and in its Invariant 2 paragraph. The roxygen of the four wrappers
      states the formula and why it never gives 0. The surrogate-testing
      vignette and the `generate_surrogate_pseudo()` roxygen state that with
      N dyads and `keep_roles = TRUE` the smallest single-dyad pseudo-dyad
      p-value is 1 / N, so p < .05 needs at least 21 dyads. The wcc and wdtw workflow vignettes no
      longer report a p-value as `< .001` or `< 0.001`. NEWS.md has an
      entry. README.md is re-knit from README.Rmd.
- [ ] AC6: `devtools::document()` leaves no diff, `devtools::test()`
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
- [ ] T4: As a finding aid, run `grep -rniE 'p-value|p_value|p value|tail'`
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

## Decisions
<!-- owner: implement / review · append-only; milestone-local. -->

## Review
<!-- owner: review · exclusive -->

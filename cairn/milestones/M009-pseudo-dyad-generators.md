<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. The one size check that can fail is
     cairn_validate's <150 over the plan-owned body. -->
# M009: Pseudo-dyad surrogate generators

- **Status:** in-progress   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate; M<xx>, M<yy> or — -->
- **Driving RR:** —   <!-- owner: plan · create/amend-via-gate -->
- **Principles touched:** —   <!-- owner: plan · create/amend-via-gate; bsync's numbered invariants live in CLAUDE.md (2, 3, 6, 7 bind here), no DESIGN IP/GP numbering yet -->
- **Resolves:** —   <!-- owner: plan · create/amend-via-gate -->
- **Surface tier:** user-facing — two new exported generator functions   <!-- owner: plan · create/amend-via-gate -->
- **Branch/PR:** m009-pseudo-dyad-generators   <!-- owner: implement (branch) / review (PR URL) · create -->

## Goal
<!-- owner: plan · create; a wrong goal returns to plan, never edited in place -->

Add two exported generators that build pseudo-dyads (one partner's series
paired with a partner from a different dyad), so users can test synchrony
against the between-dyad pseudo-synchrony null.

## Scope
<!-- owner: plan · create/amend-via-gate -->

**In:** `generate_surrogate_pseudo()`, a per-dyad partner matrix that the
four `*_surrogate()` wrappers accept. `generate_pseudo_dyads()`, the
sample-wide pseudo-dyad list for comparing real dyads to pseudo-dyads. Both
read the `dyad_list` form that `autotune_wcc()` already reads (via
`.extract_xy()`, `R/autotune.R:344`). Both crop start-aligned and keep
roles by default (`keep_roles = TRUE`). Source note for the rMEA `shuffle()`
convention. Docs, pkgdown, NEWS, and a surrogate-testing vignette section.

**Out:** pseudo-dyads as a `surrogate_method` in `synchrony_multiverse()`
and `autotune_wcc()` go to a new candidate row (this commit). IAAFT and
segment shuffling stay on the "Expanded surrogate generators" candidate row
(narrowed in this commit). A sample-level test function that compares real
dyads to pseudo-dyads stays with the "Group-level / multivariate workflow"
candidate row. This milestone adds no multi-dyad example dataset. Examples
build dyad lists from `sim_dyad` columns.

## Acceptance criteria
<!-- owner: plan · create/amend-via-gate; review reads, never reinterprets. -->

- [x] AC1: `generate_surrogate_pseudo(dyad_list, dyad, n_surrogates = NULL,
      keep_roles = TRUE)` is exported. It returns a numeric matrix with
      `length(y)` rows, where `y` is the series that `.extract_xy()` returns
      for `dyad_list[[dyad]]`. Every column equals the first `length(y)`
      samples of one partner series from a dyad other than the target. No
      two columns come from the same source series. With `keep_roles = TRUE`,
      every partner is the `y` of another dyad. With `FALSE`, a partner can
      be the `x` or the `y` of another dyad. If any partner is longer than
      `y`, a cli message states how many partners it cropped. A test matches
      each column against every candidate source series to identify its
      source. The probe series are pairwise distinct. The probes cover unequal
      lengths, a partner shorter than the target, data-frame elements, named
      and unnamed list elements, N = 2 and an odd N, both `keep_roles`
      values, and both `n_surrogates = NULL` and an integer.
- [x] AC2: `generate_surrogate_pseudo()` excludes partner series shorter
      than the target's `y`. When it excludes at least one, a cli warning
      states how many. When it excludes none, no warning fires. With
      `n_surrogates = NULL`, the column count equals the eligible partner
      count. With an integer `n_surrogates`, it samples that many partners
      without replacement. The no-eligible-partner check runs before the
      `n_surrogates` count check. Tests assert the warning and its count,
      the `NULL` column count, and fire each abort branch with a message
      matcher: fewer than two dyads, `dyad` not a single valid index into
      `dyad_list`, `n_surrogates` not a single positive integer, no eligible
      partner, and `n_surrogates` greater than the eligible partner count.
- [x] AC3: `generate_pseudo_dyads(dyad_list, n_pairs = NULL,
      keep_roles = TRUE)` is exported. It returns a list of pseudo-dyads.
      Each one is a list with named elements `x` and `y`. Its two series come
      from different dyads. Each series is the first `min(length)` samples
      of its source. If any pair is cropped, a cli message states how many.
      The attribute `"sources"` is a data frame with columns `x_dyad`,
      `x_role`, `y_dyad`, and `y_role`, one row per pseudo-dyad. With
      `n_pairs = NULL`, it returns every pairing. For `keep_roles = TRUE`
      that is the N(N-1) ordered pairs (`x` of dyad i, `y` of dyad j,
      i != j). For `keep_roles = FALSE` that is the 2N(N-1) unordered pairs
      of series from different dyads. Tests assert that the set of recorded
      source rows equals the enumerated set with no duplicates, for N = 2,
      3, and 4 under both `keep_roles` values. Tests assert that the samples
      of each element match its recorded sources. The probes cover unequal
      lengths and data-frame and list elements. With an integer `n_pairs`,
      tests assert the length equals `n_pairs` and the source rows are
      distinct. Tests fire each abort branch with a message matcher: fewer
      than two dyads, `n_pairs` not a single positive integer, and `n_pairs`
      greater than the pairing count.
- [x] AC4: Both generators draw from the session RNG and do not reseed it.
      For each generator, a test shows that two calls under the same
      `set.seed()` return identical output when they subsample. A test also
      shows that two consecutive calls without reseeding return different
      output when they subsample.
- [x] AC5: A test builds a dyad list of equal-length, NA-free series from
      `sim_dyad` columns. It passes one `generate_surrogate_pseudo()` matrix
      to each of `wcc_surrogate()`, `wdtw_surrogate()`,
      `wgranger_surrogate()`, and `wphase_surrogate()`, with `window_size`
      and `lag_max` valid for that length. It asserts that `p_value` (both
      `p_value_xy` and `p_value_yx` for Granger) is in [0, 1]. It asserts that
      `n_surrogates` equals the column count of the matrix.
- [x] AC6: The roxygen of both functions states the rationale for each
      default, including why `keep_roles` defaults differ from rMEA. It has
      runnable `@examples` and `@seealso` links to the surrogate wrappers and
      the existing generators. Both functions are in the surrogate section of
      `_pkgdown.yml`. NEWS.md has an entry. The surrogate-testing vignette
      gets a pseudo-dyad section with a worked example. The section states
      which null hypothesis the pseudo-dyad null tests, and cites the
      published convention source.
- [x] AC7: `devtools::document()` produces no diff against the committed
      tree. `devtools::test()` passes. `devtools::check()` reports 0 errors
      and 0 warnings, and the milestone file triages any NOTE.

## Coverage
<!-- owner: plan · create/amend-via-gate -->

- AC1 → T2, T3
- AC2 → T2, T3
- AC3 → T4
- AC4 → T2, T3, T4
- AC5 → T5
- AC6 → T1, T6, T8
- AC7 → T7

## Tasks
<!-- owner: plan (create) / implement (check-off, minor edits) -->

- [x] T1: Write the source note for the rMEA `shuffle()` convention
      (CRAN source `R/rMEA_rand.R`, cited as Kleinbub & Ramseyer 2020) in
      `cairn/references/` with its `INDEX.md` line. Record the pairing,
      role, sampling, and cropping rules, and where bsync departs from them.
- [x] T2: Tests first for `generate_surrogate_pseudo()` in
      `tests/testthat/test-surrogate-pseudo.R` (AC1, AC2, AC4). Use pairwise
      distinct seeded series. Assert cli output with `expect_message()` or
      `expect_warning()` (LESSONS 2026-08-29).
- [x] T3: Implement `generate_surrogate_pseudo()` in
      `R/surrogate_generation.R`. Reuse `.extract_xy()` and move it to a
      shared helper file if both new functions need it.
- [x] T4: Tests first, then implement `generate_pseudo_dyads()` (AC3, AC4).
- [x] T5: Wrapper integration test over all four `*_surrogate()` functions
      (AC5).
- [x] T6: Docs. Roxygen for both functions, `_pkgdown.yml` surrogate
      section, NEWS.md entry, surrogate-testing vignette section,
      `inst/WORDLIST`. Update `cairn/DESIGN.md` §6 generator list and the
      CLAUDE.md out-of-scope line so they say the pseudo-dyad generators
      shipped.
- [x] T7: Gate. Run `devtools::document()`, `devtools::test()`,
      `devtools::check()`, `air format`, and `lintr::lint_package()`.
- [x] T8: Review return 1. Add @seealso links from `generate_pseudo_dyads()`
      to the four surrogate wrappers and the existing generators (AC6).

## Work log
<!-- owner: any skill · append-only; one line per entry; absolute dates. -->

- 2026-10-05: created by /milestone-plan.
- 2026-10-05: collision sweep: absorbed the pseudo-dyad part of the "Expanded surrogate generators" candidate row. IAAFT and segment shuffling stay on that row (narrowed). No DECISIONS entries, no archive overlap (M008 is wphase), 0 open issues, 0 open PRs.
- 2026-10-05: no `Depends on:` the multi-dyad workflow row, although the candidate row and CLAUDE.md said so, because `autotune_wcc()` already defines the `dyad_list` input these generators read.
- 2026-10-05: criteria audit (full mode, fresh Opus reader) returned 14 findings across AC1–AC7. All were accepted and fixed in the written wording: `.extract_xy()` definition of `y`, cli message on cropping (Invariant 3), wider probe axes and distinct probe series, warning only when excluding, abort-check order, `"sources"` attribute shape, set-equality for pairings, RNG test replaces a `set.seed` grep, Granger's two p-values, a named `sim_dyad` fixture, `keep_roles` rationale, NOTE triage.
- 2026-10-05: question set: functions: both generators. Alignment: start-aligned crop. Roles: keep by default, with a `keep_roles` option. Multiverse integration: later, as a candidate row.
- 2026-10-05: plan gate chose start-aligned cropping over random-offset cropping because it keeps task-locked timing in the null and matches rMEA; falsified by evidence that start-aligned pseudo-dyads keep within-dyad coupling or give miscalibrated p-values on uncoupled simulated dyads.
- 2026-10-05: plan gate chose `keep_roles = TRUE` as the default over rMEA's role-mixing default because bsync's lead-lag sign depends on which series is `x`; falsified by evidence that users mostly analyze dyads with interchangeable partners.
- 2026-10-05: plan chose a `"sources"` data-frame attribute over per-element source fields because it keeps each element in the plain `list(x, y)` dyad form; falsified by a downstream need to carry sources through `lapply()`.
- 2026-10-05: plan chose to abort when `n_surrogates` exceeds the eligible count, not to sample with replacement as the circular generator does, because duplicate partners add no null information.
- 2026-10-05: no new example dataset: examples and the AC5 fixture build dyad lists from `sim_dyad` columns. The vignette simulates a small dyad list inline.
- 2026-10-05: branch m009-pseudo-dyad-generators cut from main at 39ca324; status in-progress.
- 2026-10-05: T1 done: cairn/references/kleinbub2020.md from the rMEA 1.2.2 source (shuffle(), unequalCbind()). It confirms the start-aligned crop to the shorter series and the N(N-1) and 2N(N-1) pair counts. Departures recorded there: keep_roles default, no na.omit().
- 2026-10-05: T2+T3 done: generate_surrogate_pseudo() in R/surrogate_generation.R with tests in test-surrogate-pseudo.R (11 tests). Full suite: 0 failures. Planted defects went red: end-aligned crop (3 tests fail), target's own series in the pool (8 tests fail).
- 2026-10-05: implementation choices: .extract_xy() stays in R/autotune.R (package-internal, reachable from both files). Added a numeric-series abort and a keep_roles TRUE/FALSE abort beyond AC2's list. The matrix also carries a "sources" attribute (dyad, role) matching AC3's form. NULL n_surrogates keeps pool order and draws no random numbers.
- 2026-10-05: T4 done: generate_pseudo_dyads() with tests that compare the recorded pairings to an independent enumeration for N = 2, 3, 4 under both keep_roles values. Full suite: 1512 expectations, 0 failures. Planted defects went red: same-dyad pairs with roles kept (4 tests fail), same-dyad pairs with roles mixed (3), end-aligned crop (1).
- 2026-10-05: T5 done: wrapper test on a 3-dyad list from sim_dyad's axes (first 400 samples, window_increment 20 to keep WDTW fast). All four wrappers return p-values in [0, 1] and n_surrogates = 2. A further test runs wcc() on a generate_pseudo_dyads() element.
- 2026-10-05: T6 done: roxygen, _pkgdown.yml, NEWS.md, WORDLIST, DESIGN §2 and §6, the CLAUDE.md out-of-scope line, and the DESCRIPTION Description field. The vignette gets section 1.3 and a worked example (section 3). Its rendered output shows 9 of 10 task-only dyads significant against circular shifts and 2 of 10 against pseudo-dyads, and real dyads scoring higher than pseudo-dyads. The prose was written against that render.
- 2026-10-05: found that roxygen markdown is off package-wide, so existing `[fn()]` links render as raw text. Used a per-block `@md` on the two new functions to make their @seealso links render (AC6). Package-wide fix recorded as a candidate row, not done here.
- 2026-10-05: T7 gate (first pass, pre-claim-audit code): air format applied to the two new files. lintr: no hits in new code (the object_usage hits are package-wide false positives, bsync not installed). document() no diff except man/bsync-package.Rd from the DESCRIPTION edit. check_pkgdown(): no problems. check(): 0 errors, 0 warnings, 0 notes.
- 2026-10-05: claim audit: 66 claims read, 3 corrected — R/surrogate_generation.R, vignettes/surrogate-testing.Rmd. Crop message counted partners before the n_surrogates draw (code fixed). Vignette implied both functions crop to the shorter series (wording fixed). keep_roles = FALSE gave x = 1y for a (1y, 2x) pair where rMEA gives x = 2x (pool reordered to rMEA's 1x..Nx, 1y..Ny, roxygen states it). Two regression tests added; both fail on the pre-fix code.
- 2026-10-05: T7 done on the final code (301e560): check() 0 errors, 0 warnings, 0 notes; document() no diff. Status set to review.
- 2026-10-05: review return 1: AC6 fails as written. The `generate_pseudo_dyads()` roxygen @seealso has no links to the surrogate wrappers (`wcc_surrogate()`, `wdtw_surrogate()`, `wgranger_surrogate()`, `wphase_surrogate()`) or the existing generators (`generate_surrogate_circular()`, `generate_surrogate_phase()`). AC1–AC5 and AC7 passed with evidence. Status back to in-progress.
- 2026-10-05: T8 added for review return 1 (minor amendment; Coverage AC6 → T1, T6, T8). Done: generate_pseudo_dyads() @seealso now links the four `*_surrogate()` wrappers and both existing generators (12 links in the Rd). Its one new claim, that the wrappers accept the per-dyad matrix, is the behavior test 17 asserts. Full suite: 1529 expectations, 0 failures. Status set to review.
- 2026-10-05: review return 2: consistency gate FAIL, `cairn_validate` weight caps. The plan-owned body is 150 lines (cap <150) after T8 was added. Heaviest sections: Acceptance criteria 71, Tasks 27. Status back to in-progress.

## Decisions
<!-- owner: implement / review · append-only; milestone-local -->

## Review
<!-- owner: review · exclusive -->

Pass 1 (2026-10-05, branch head bb1adf9). `testthat::test_file("tests/testthat/test-surrogate-pseudo.R")`: 20 tests, all passing.

- AC1 evidence: `export(generate_surrogate_pseudo)` in NAMESPACE. Source identity is matched against every candidate series (`column_sources()`, pairwise distinct `rnorm` probes) in tests 1 (N = 5, df/named/unnamed/reversed-named elements, unequal lengths, one shorter partner, keep_roles TRUE), 2 (same, keep_roles FALSE), 3 (y read by name), 4 (N = 2), 5 (integer n_surrogates), 7. The crop message is asserted in tests 1, 2, and 19 (count of returned partners). Pass.
- AC2 evidence: excluded-count warning asserted in tests 1 and 2. No warning when nothing is excluded: tests 4 and 8. NULL column count equals the eligible count: tests 1, 2, 6. Integer draw gives distinct partners: test 5. Abort branches with message matchers: test 9 (fewer than two dyads, invalid `dyad`, invalid `n_surrogates`, count over eligible) and test 10 (no eligible partner, which fires first when `n_surrogates` also exceeds zero). Pass.
- AC3 evidence: `export(generate_pseudo_dyads)` in NAMESPACE. Test 12 compares the recorded source set to an independent enumeration with no duplicates, for N = 2, 3, 4 under both keep_roles values (N(N-1) and 2N(N-1)). It also checks the `"sources"` column names, that x_dyad != y_dyad, the `x`/`y` element names, and that samples match their sources (unequal lengths; df and list elements). Test 13: crop message. Test 14: integer `n_pairs` length and distinct rows. Test 15: the three abort branches. Pass.
- AC4 evidence: same seed gives identical output, and consecutive unseeded calls differ, while subsampling: test 11 (generate_surrogate_pseudo) and test 16 (generate_pseudo_dyads). Pass.
- AC5 evidence: test 17 builds a 3-dyad list from sim_dyad x/y/z columns (400 samples, NA-free) and runs all four wrappers. `p_value` is in [0, 1] (both Granger p-values) and `n_surrogates` = 2 = ncol. Pass.
- AC6 evidence: FAIL. `man/generate_pseudo_dyads.Rd` links to autotune_wcc, generate_surrogate_pseudo, wcc, wdtw, wgranger, and wphase. It does not link to the surrogate wrappers (`*_surrogate()`) or the existing generators (`generate_surrogate_circular()`, `generate_surrogate_phase()`), as the criterion requires for both functions. The other items pass: the rationale for every default including the keep_roles departure from rMEA (both Rd files), `_pkgdown.yml` lines 62–63, NEWS.md lines 22 and 30, and vignette section 1.3 "Which null it tests" with the Kleinbub & Ramseyer (2020) citation and the section 3 worked example.
- AC7 evidence: at bb1adf9, `devtools::document()` leaves man/ and NAMESPACE unchanged. `devtools::check()` reports 0 errors, 0 warnings, 0 notes, and it runs the full test suite, the examples, and the vignette. Pass.
- Consistency gate: `cairn_validate.py` all checks passed. `pkgdown::check_pkgdown()` no problems. README.Rmd and README.md are unchanged on the branch. No new top-level files. The NEWS.md entry is present with no milestone numbers. Principle change: none, so `cairn_impact` is skipped.

Pass 2 (2026-10-05, branch head 6247c38, after review return 1). `test-surrogate-pseudo.R`: 20 tests, 357 expectations, all passing. The only code change since pass 1 is roxygen @seealso text, so the AC1–AC5 evidence stands, re-run at this head.
- AC6 evidence (pass 2): `man/generate_pseudo_dyads.Rd` now links wcc_surrogate, wdtw_surrogate, wgranger_surrogate, wphase_surrogate, generate_surrogate_circular, and generate_surrogate_phase (also autotune_wcc, generate_surrogate_pseudo, wcc, wdtw, wgranger, wphase). `man/generate_surrogate_pseudo.Rd` links the same wrappers and generators. The other AC6 items are unchanged from pass 1. Pass.
- AC7 evidence (pass 2): at 6247c38, `devtools::check()` reports 0 errors, 0 warnings, 0 notes. `document()` gives no diff.
- spawned: diff-bug, blame-history, prior-review
- diff-bug #1: with N dyads, a single-dyad pseudo-dyad test gives p = 0 with chance about 1/N under the null (p = b/n, no +1), and the vignette did not say so — fix now (roxygen "How many surrogates" and vignette section 3 state the 1/N rate and point small samples to the sample-level comparison), fixed 7302851. The wrapper formula change goes to the candidate row "Add-one surrogate p-value".
- diff-bug #2: a zero-length target `y` returned a 0-row matrix instead of aborting — fix now (abort in `.pseudo_extract_dyads()`, test added), fixed 7302851.
- diff-bug #3: `.extract_xy()` partial-matches names and silently falls back to position — follow-up, a pre-existing helper shared with `autotune_wcc()`: candidate row "Strict dyad extraction in `.extract_xy()`".
- diff-bug #4: no test of NA pass-through or integer input — fix now (test added), fixed 7302851.
- diff-bug #5: no test of `n_surrogates` / `n_pairs` equal to the full count — fix now (test added), fixed 7302851.
- diff-bug #6: the `"sources"` attribute is dropped on subsetting — fix now (documented in both @return), fixed 7302851.
- diff-bug #7: `dyad_list` names are not carried into `sources` — follow-up: candidate row "Carry `dyad_list` names into the pseudo-dyad `"sources"` attribute".
- diff-bug #8: the docs did not state the shared sampling-rate and shared-onset assumption behind start alignment — fix now (both roxygen blocks and the vignette), fixed 7302851.
- diff-bug #9: a bare data frame as `dyad_list` got a misleading "at least two dyads" error — fix now (separate message, test added), fixed 7302851.
- diff-bug #10: `.extract_xy()` errors do not name the dyad index — follow-up: same candidate row as #3.
- diff-bug #11: `generate_surrogate_pseudo()` does not check the target's `x` length against `y` — reject (false as a defect): the generator's contract covers `y` only, and every wrapper aborts with "x and y must be the same length" before computing.
- diff-bug #12: the `sources`-to-columns check ran only for keep_roles = FALSE — fix now (test added for TRUE), fixed 7302851.
- diff-bug #13: DESCRIPTION labels all generators "(pseudo-synchrony)" and the line is long — reject (false): the "(pseudo-synchrony)" label predates M009, and DESIGN §1 uses the term for all surrogate nulls. Line length is style.
- diff-bug #14: DESIGN §6 wrapper list omits `wphase_surrogate()` — fix now (pre-existing, corrected in place), fixed 7302851.
- blame-history #1: README vignette row and surrogate references still described circular-shift and phase only — fix now (README.Rmd row and rMEA reference, README.md re-knit), fixed 7302851.
- blame-history #2: DESIGN §14 #8 still said "roadmap M12" — fix now (corrected in place), fixed 7302851.
- blame-history #3: DESIGN §6 lead sentence said every generator preserves each series' own structure — fix now (qualified for pseudo-dyads), fixed 7302851.
- blame-history #4: the vignette's "Decouple" step and "Interpretation" bullet assumed within-dyad nulls — fix now (both qualified), fixed 7302851.
- blame-history #5: CLAUDE.md reverses the old "depends on M11's dyad_list" claim — reject (false as a defect): the reversal matches the code (`R/autotune.R:63`) and no D-entry records the old claim.
- blame-history #6: DESCRIPTION edit consistent with man/bsync-package.Rd — noted, no issue.
- blame-history #7: NEWS.md additions alter no existing text — noted, no issue.
- prior-review #1: the `generate_pseudo_dyads()` example runs `wcc()` on six pseudo-dyads outside `\donttest{}` — reject (false): measured at 0.05 s elapsed.
- After the fixes: `test-surrogate-pseudo.R` passes 24 tests. The full suite passes 1561 expectations with 0 failures. Spell check is clean. `document()` gives no diff.

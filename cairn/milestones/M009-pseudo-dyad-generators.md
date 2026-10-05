<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. The one size check that can fail is
     cairn_validate's <150 over the plan-owned body. -->
# M009: Pseudo-dyad surrogate generators

- **Status:** planned   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate; M<xx>, M<yy> or — -->
- **Driving RR:** —   <!-- owner: plan · create/amend-via-gate -->
- **Principles touched:** —   <!-- owner: plan · create/amend-via-gate; bsync's numbered invariants live in CLAUDE.md (2, 3, 6, 7 bind here), no DESIGN IP/GP numbering yet -->
- **Resolves:** —   <!-- owner: plan · create/amend-via-gate -->
- **Surface tier:** user-facing — two new exported generator functions   <!-- owner: plan · create/amend-via-gate -->
- **Branch/PR:** —   <!-- owner: implement (branch) / review (PR URL) · create -->

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

- [ ] AC1: `generate_surrogate_pseudo(dyad_list, dyad, n_surrogates = NULL,
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
- [ ] AC2: `generate_surrogate_pseudo()` excludes partner series shorter
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
- [ ] AC3: `generate_pseudo_dyads(dyad_list, n_pairs = NULL,
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
- [ ] AC4: Both generators draw from the session RNG and do not reseed it.
      For each generator, a test shows that two calls under the same
      `set.seed()` return identical output when they subsample. A test also
      shows that two consecutive calls without reseeding return different
      output when they subsample.
- [ ] AC5: A test builds a dyad list of equal-length, NA-free series from
      `sim_dyad` columns. It passes one `generate_surrogate_pseudo()` matrix
      to each of `wcc_surrogate()`, `wdtw_surrogate()`,
      `wgranger_surrogate()`, and `wphase_surrogate()`, with `window_size`
      and `lag_max` valid for that length. It asserts that `p_value` (both
      `p_value_xy` and `p_value_yx` for Granger) is in [0, 1]. It asserts that
      `n_surrogates` equals the column count of the matrix.
- [ ] AC6: The roxygen of both functions states the rationale for each
      default, including why `keep_roles` defaults differ from rMEA. It has
      runnable `@examples` and `@seealso` links to the surrogate wrappers and
      the existing generators. Both functions are in the surrogate section of
      `_pkgdown.yml`. NEWS.md has an entry. The surrogate-testing vignette
      gets a pseudo-dyad section with a worked example. The section states
      which null hypothesis the pseudo-dyad null tests, and cites the
      published convention source.
- [ ] AC7: `devtools::document()` produces no diff against the committed
      tree. `devtools::test()` passes. `devtools::check()` reports 0 errors
      and 0 warnings, and the milestone file triages any NOTE.

## Coverage
<!-- owner: plan · create/amend-via-gate -->

- AC1 → T2, T3
- AC2 → T2, T3
- AC3 → T4
- AC4 → T2, T3, T4
- AC5 → T5
- AC6 → T1, T6
- AC7 → T7

## Tasks
<!-- owner: plan (create) / implement (check-off, minor edits) -->

- [ ] T1: Write the source note for the rMEA `shuffle()` convention
      (CRAN source `R/rMEA_rand.R`, cited as Kleinbub & Ramseyer 2020) in
      `cairn/references/` with its `INDEX.md` line. Record the pairing,
      role, sampling, and cropping rules, and where bsync departs from them.
- [ ] T2: Tests first for `generate_surrogate_pseudo()` in
      `tests/testthat/test-surrogate-pseudo.R` (AC1, AC2, AC4). Use pairwise
      distinct seeded series. Assert cli output with `expect_message()` or
      `expect_warning()` (LESSONS 2026-08-29).
- [ ] T3: Implement `generate_surrogate_pseudo()` in
      `R/surrogate_generation.R`. Reuse `.extract_xy()` and move it to a
      shared helper file if both new functions need it.
- [ ] T4: Tests first, then implement `generate_pseudo_dyads()` (AC3, AC4).
- [ ] T5: Wrapper integration test over all four `*_surrogate()` functions
      (AC5).
- [ ] T6: Docs. Roxygen for both functions, `_pkgdown.yml` surrogate
      section, NEWS.md entry, surrogate-testing vignette section,
      `inst/WORDLIST`. Update `cairn/DESIGN.md` §6 generator list and the
      CLAUDE.md out-of-scope line so they say the pseudo-dyad generators
      shipped.
- [ ] T7: Gate. Run `devtools::document()`, `devtools::test()`,
      `devtools::check()`, `air format`, and `lintr::lint_package()`.

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

## Decisions
<!-- owner: implement / review · append-only; milestone-local -->

## Review
<!-- owner: review · exclusive -->

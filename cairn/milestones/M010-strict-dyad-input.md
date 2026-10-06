<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. The one size check that can fail is
     cairn_validate's <150 over the plan-owned body. -->
# M010: Strict dyad input and named pseudo-dyad sources

- **Status:** review   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate; M<xx>, M<yy> or — -->
- **Driving RR:** —   <!-- owner: plan · create/amend-via-gate -->
- **Principles touched:** —   <!-- owner: plan · create/amend-via-gate; bsync's numbered invariants live in CLAUDE.md (3 binds here), no DESIGN IP/GP numbering yet -->
- **Resolves:** —   <!-- owner: plan · create/amend-via-gate -->
- **Surface tier:** user-facing — changes how exported functions read `dyad_list` and what the pseudo-dyad generators return   <!-- owner: plan · create/amend-via-gate -->
- **Branch/PR:** m010-strict-dyad-input   <!-- owner: implement (branch) / review (PR URL) · create -->

## Goal
<!-- owner: plan · create; a wrong goal returns to plan, never edited in place -->

Make every function that reads a `dyad_list` identify each dyad's series by
exact name or by position, never by partial match or a silent role swap,
with errors that name the dyad, and carry `dyad_list` names into the
pseudo-dyad generators' output.

## Scope
<!-- owner: plan · create/amend-via-gate -->

**In:** A strict `.extract_xy()` (`R/autotune.R:344`), shared by
`autotune_wcc()` (`R/autotune.R:177`) and `.pseudo_extract_dyads()`
(`R/surrogate_generation.R:474`). One rule for list and data-frame dyads:
read by the exact names `x` and `y` when both occur once, read the first two
elements when neither occurs, else abort. Every extraction abort names the
dyad's index in `dyad_list`. `autotune_wcc()` reads every dyad before it
samples `n_tune_dyads`. The `"sources"` attribute of
`generate_surrogate_pseudo()` and `generate_pseudo_dyads()` gets name
columns. If `dyad_list` is named, the matrix columns and list elements also
get names. Roxygen and NEWS.

**Out:** Pseudo-dyads as a `surrogate_method` in `synchrony_multiverse()`
stay on their candidate row. The add-one p-value is M011 (planned in the
same commit). Dyad names in `autotune_wcc()` output are not part of either
candidate row and are not added here.

## Acceptance criteria
<!-- owner: plan · create/amend-via-gate; review reads, never reinterprets. -->

- [ ] AC1: For a dyad that is a list or a data frame, `.extract_xy()`
      applies one rule. If the element (column) names contain `x` exactly
      once and `y` exactly once, it returns those two elements. If the names
      contain neither `x` nor `y`, or there are no names, it returns the
      first two elements, and a dyad with fewer than two elements is an
      abort. Any other count of exact `x` or `y` names is an abort. No name
      is matched partially. Tests through `generate_surrogate_pseudo()`
      assert which source series each returned column holds for: `list(y =
      b, x = a)`, a data frame with columns `time`, `x`, `y`, `list(yy = b,
      xx = a)` (read by position, so `x` is `b`), and an unnamed list. Tests
      fire the abort branches with a message matcher for: `list(y = a, a2 =
      b)`, `list(x = a, b2 = b)`, `list(x = a, x = b, y = c)`, a data frame
      with two columns named `x` (`check.names = FALSE`), a one-column data
      frame, and a non-list dyad.
- [ ] AC2: Every abort raised by `.extract_xy()` names the dyad's index in
      the user's `dyad_list`. `autotune_wcc()` reads every dyad before it
      samples. Thus a dyad that `.extract_xy()` rejects aborts the call also
      when the sample leaves it out. Tests assert the index in the message
      for a rejected dyad at position 3 of 4 through
      `generate_surrogate_pseudo()`, through `generate_pseudo_dyads()`, and
      through `autotune_wcc()` called with `n_tune_dyads = 2`. Each path is
      probed with two forms: `list(y = a, a2 = b)` and a one-column data
      frame.
- [ ] AC3: A `dyad_list` counts as named when `names(dyad_list)` is
      non-NULL and every name is non-empty, not `NA`, and unique. With a
      named list, the `"sources"` data frame of `generate_surrogate_pseudo()`
      gains a character column `name`, the source dyad's name, and the
      matrix column names are `"<name>:<role>"`. The `"sources"` data frame
      of `generate_pseudo_dyads()` gains `x_name` and `y_name`, and each list
      element is named `"<x_name>:<x_role>|<y_name>:<y_role>"`. With
      `names(dyad_list)` NULL, the name columns are `NA_character_` and the
      matrix and list carry no names. Any other naming (some names empty or
      `NA`, or a name repeated) is an abort in both generators. The integer
      columns `dyad`, `x_dyad`, and `y_dyad` keep their M009 meaning. Tests
      cover a named and an unnamed list for both generators under both
      `keep_roles` values and under an integer `n_surrogates` / `n_pairs`,
      and assert that each name matches the row's integer index. Tests fire
      the empty-name, `NA`-name, and repeated-name aborts with a message
      matcher. Names that contain `:` or `|` are accepted as given.
- [ ] AC4: The `dyad_list` roxygen of `generate_surrogate_pseudo()` (which
      `generate_pseudo_dyads()` inherits) and of `autotune_wcc()` states the
      AC1 rule. The `@return` of both generators states the name columns and
      when names appear. NEWS.md has an entry for the stricter reading and the
      names. The entry names the input forms whose result changes without
      an error: a data frame or list with elements named `x` and `y` that
      are not its first two in that order (for example columns `time`, `x`,
      `y`), and a list whose names only start with `x` or `y` (for example
      `list(yy =, xx =)`), which is now read by position. `devtools::document()` leaves no diff, `devtools::test()` passes,
      and `devtools::check()` reports 0 errors and 0 warnings.

## Coverage
<!-- owner: plan · create/amend-via-gate -->

- AC1 → T1
- AC2 → T1
- AC3 → T2, T4
- AC4 → T3

## Tasks
<!-- owner: plan (create) / implement (check-off, minor edits). -->

- [x] T1: Tests first in `tests/testthat/test-surrogate-pseudo.R` and
      `test-autotune.R`. Rewrite `.extract_xy(dyad, index)` to the AC1 rule
      with `[[` on exact names, never `$`. Pass the index from
      `.pseudo_extract_dyads()`. In `autotune_wcc()`, extract every dyad
      before the `n_tune_dyads` sample (`R/autotune.R:161-177`) and run the
      multiverse on the extracted series.
- [x] T2: Tests first. Add a name check to `.pseudo_extract_dyads()` (or a
      sibling helper) that returns the names or `NULL`. Add the name columns
      to both `"sources"` data frames after any subsample, and set the matrix
      column names and list element names.
- [x] T3: Update the roxygen (`@md` is already on both generators; add it to
      `autotune_wcc()` only if the block gains markdown), run
      `devtools::document()`, add the NEWS entry, run `devtools::test()` and
      `devtools::check()`, style with `air format`.
- [x] T4: (review return 1) Extend the "unnamed dyad_list" test in
      `test-surrogate-pseudo.R` to an integer `n_surrogates` and `n_pairs`
      draw for both generators under both `keep_roles` values.

## Work log
<!-- owner: any skill · append-only; one line per entry; absolute dates. -->

- 2026-10-05: created by /milestone-plan.
- 2026-10-05: question set: what to work on — "Pseudo-dyad cleanup". The add-one p-value answer ("Yes, add to this plan") went to M011.
- 2026-10-05: lineage: absorbs the candidate rows "Strict dyad extraction in `.extract_xy()`" and "Carry `dyad_list` names into the pseudo-dyad `sources` attribute" (both M009 review); both rows removed in the plan commit.
- 2026-10-05: inbox sweep: 0 open issues, 0 open PRs.
- 2026-10-05: split: M010 covers dyad reading and names; the add-one p-value is M011 (planned now, no dependency, ROADMAP order puts M010 first).
- 2026-10-05: plan gate chose exact-names-else-position over aborting on any naming other than `x`/`y` because the abort breaks the documented positional form for named lists such as `list(p1 =, p2 =)`; falsified by a user report of a named list read by position against intent.
- 2026-10-05: plan gate chose to read every dyad before sampling in `autotune_wcc()` over mapping sampled positions back to indices because it makes the index exact and aborts loudly on any malformed dyad; falsified by a `dyad_list` large enough that up-front extraction costs noticeable time or memory.
- 2026-10-05: plan gate chose added `name` columns over converting the integer `dyad` columns to character because it keeps the M009 columns unchanged; falsified by users needing the name as the only key in `"sources"`.
- 2026-10-05: criteria audit (full mode, fresh Opus reader) returned 4 findings on M010. Fixed: AC1 partial-match probe swapped to `list(yy = b, xx = a)`, which tells the two readings apart, and abort probes widened (`list(x = a, b2 = b)`, duplicate-`x` data frame). Fixed: AC2 narrowed to dyads that `.extract_xy()` rejects, two forms per path. Fixed: AC3 gains the `NA`-name probe, and names with `:` or `|` are accepted. Fixed: AC4 NEWS must name the inputs whose result changes without an error.
- 2026-10-05: plan gate chose a NEWS entry over a per-call cli notice for inputs now read by name instead of by position, because after the upgrade the notice repeats on every correct call; falsified by users who miss the change and report shifted results.
- 2026-10-05: implement started on branch `m010-strict-dyad-input`, cut from main at 3115dba.
- 2026-10-05: T1 done. `.extract_xy(dyad, index)` applies the exact-name-else-position rule with abort messages that name the dyad index; `autotune_wcc()` extracts every dyad before sampling. New tests saw red first (7 failures), then the full suite passed (1579 tests). `air` restyled untouched code in `R/autotune.R` and `test-autotune.R`, which were never air-clean, so that restyle was dropped to keep the diff focused.
- 2026-10-05: T2 done. New helper `.pseudo_dyad_names()` returns the names or NULL and aborts on empty, NA, or repeated names. Both generators add the name columns after any draw. On the pre-T2 code the new tests failed (45 failures), and on the new code the full suite passed (1710 tests). The M009 test that pinned the 4-column `generate_pseudo_dyads()` "sources" schema now pins the 6-column schema (AC3).
- 2026-10-05: implementation choice: in `generate_pseudo_dyads()` the "sources" data frame is rebuilt with the name columns placed after each role column, which gives the AC3 column order.
- 2026-10-05: T3 done. Roxygen for `dyad_list` (both generators and `autotune_wcc()`) and the two `@return` blocks updated. `autotune_wcc()` stays without `@md` and uses `\code{}`. NEWS has a "Stricter `dyad_list` reading" section. `devtools::document()` is clean on a second run, `devtools::check()` gives 0 errors, 0 warnings, 0 notes, and `spelling::spell_check_package()` finds nothing. `R/autotune.R` is still not air-clean, as it was on main.
- 2026-10-05: the T2 commit swept in testthat failure files from the red run (`tests/testthat/_problems/`, `testthat-problems.rds`). They are removed and added to `.gitignore`.
- 2026-10-05: claim audit: 44 claims read, 4 corrected — NEWS.md, R/autotune.R, R/surrogate_generation.R, man/*.Rd (the silent-change bullet named lists where only data frames change and missed `list(xval =, other =, yval =)`; "read and checked" narrowed to names and shape; NA names added to the naming rule). The same reader re-read the corrections and found them true, plus one stale test comment, which is now fixed. The autotune test also runs with `n_tune_dyads = 4`, so it covers the case where the bad dyad is sampled.
- 2026-10-05: review return 1: AC3 requires the unnamed-list case under an integer `n_surrogates` / `n_pairs`, but the "unnamed dyad_list" test covers only the NULL draw for both generators. AC1, AC2, and AC4 probes are all present, and the two test files pass (42 tests, 0 failed).
- 2026-10-05: minor amendment: T4 added for review return 1, Coverage AC3 → T2, T4.
- 2026-10-05: T4 done. The unnamed-list test now loops over a NULL and an integer draw (`n_surrogates = 1`, `n_pairs = 3`) for both generators and both `keep_roles` values, and checks the source identity of each column or pair. The full suite passed (1806 tests). Status set to review.

## Decisions
<!-- owner: implement / review · append-only; milestone-local. -->

## Review
<!-- owner: review · exclusive -->

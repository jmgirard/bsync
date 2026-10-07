# M015: Segment-shuffling surrogate generator

- **Status:** blocked
- **Priority:** normal
- **Depends on:** M014
- **Driving RR:** —
- **Principles touched:** —
- **Resolves:** —
- **Surface tier:** user-facing — adds an exported generator and a new `surrogate_method` value
- **Branch/PR:** m015-segment-shuffle-surrogates

## Goal

Add an exported segment-shuffling surrogate generator that reorders whole segments of a series, with no segment in its original place. Offer it as a `surrogate_method` in `synchrony_multiverse()` and `autotune_wcc()`.

## Scope

**In:**
- `generate_surrogate_segment(y, segment_size, n_surrogates = 100)` in `R/surrogate_generation.R`. It cuts `y` into k = floor(length(y) / segment_size) complete segments from sample 1. This is the SUSY 0.1.0 convention (`susy()` source). Each column puts the segments in a random order in which no segment keeps its original position. The leftover tail of length(y) - k * segment_size samples stays in place at the end. The call announces the tail length.
- `segment_size` is a whole number from 2 to floor(length(y) / 2), so k is at least 2. This limit is adapted from SUSY's check, which compares seconds to rows.
- Orders are distinct, drawn without replacement. If `n_surrogates` is at least the number of possible orders, the call returns every possible order and says how many. Draws with replacement make the test reject too often when k is small.
- `NA` values move with their segments.
- `"segment"` as a `surrogate_method` value in `synchrony_multiverse()` and `autotune_wcc()`. The segment size is each cell's window size in samples. One matrix is generated per distinct window size and reused across the cells that share it. A cell whose window is above half the series length is skipped like other unusable cells, with a warning that names the window.
- Roxygen with the null this tests, the tail rule, how windows line up with segments, and the SUSY lineage. A `_pkgdown.yml` row, a vignette section, NEWS, `inst/WORDLIST`.
- Source note `cairn/references/tschacher2020.md` for the SUSY segment convention, and DESIGN §2, §6, §14 #8 updates.

**Out:**
- SUSY-exact segment pairing, a WCC-only statistic over every non-matching segment pair. The plan gate chose the reordered form. Not planned.
- IAAFT → M014.
- Pseudo-dyad surrogates in the multiverse: their existing candidate row.

## Acceptance criteria

- [ ] AC1: Tests on seeded inputs show that `generate_surrogate_segment()` returns a numeric matrix with `length(y)` rows and `n_surrogates` columns. In every column, the first k * `segment_size` samples split into k segments. Each one equals the segment of `y` at a different position. Each segment of `y` appears exactly once. The tail equals the tail of `y`. The inputs include a length that `segment_size` divides, a length that it does not divide, and a series with `NA` values. The `NA` positions move with their segments. For one case, the test writes the expected segment boundaries by hand from the `eps` loop of SUSY 0.1.0 (`R/susy.R`).
- [ ] AC2: A test file holds a plain reimplementation of the generator with explicit loops and the same sequence of random draws. The test uses at least two seeded cases, one with a tail and one without. In each case, the output equals the reference output (`expect_equal`).
- [ ] AC3: A seeded simulation test (`skip_on_cran()`) draws at least 200 independent pairs. Each partner is an AR(1) series x_n = 0.7 x_{n-1} + eta_n. The test runs `wcc_surrogate()` with at least 19 segment surrogates per pair, at least 6 segments, and `lag_max` of at least 1. Both `window_size` and `window_increment` equal `segment_size`. The proportion of pairs with p <= .05 is at most 0.09. A second run of the same design with k = 3 segments and `n_surrogates = 19` also gives a proportion of at most 0.09.
- [ ] AC4: A test shows that the tail message names the tail length. With k = 3 segments (two possible orders) and `n_surrogates = 5`, the call returns two distinct columns and its message names the number 2. Tests fire each abort branch and assert each message with a regexp. The branches are non-numeric `y`, an invalid `n_surrogates`, and a `segment_size` that is not a whole number of at least 2. The last branch is a `segment_size` above floor(length(y) / 2).
- [ ] AC5: A test runs a one-cell `synchrony_multiverse(estimator = "wcc", surrogate_method = "segment")`. Its `p` equals the `p` of `wcc_surrogate()` on `generate_surrogate_segment()` output drawn after the same `set.seed()`. That call uses a `segment_size` equal to the cell's window size, and the cell's window, lag, and increment in samples. A test shows that a grid cell whose window is above half the series length is skipped, with a warning that names that window. A test shows that `autotune_wcc(surrogate_method = "segment")` returns a result whose grid holds `"segment"`.
- [ ] AC6: The roxygen of `generate_surrogate_segment()` states the null it tests, the tail rule, and the SUSY lineage. It states the condition under which windows line up with segments. It cites Tschacher & Meier (2020) with its DOI and has a runnable `sim_dyad` example. The `surrogate_method` docs of both functions name `"segment"` and state that the segment size is the cell's window size. The surrogate-testing vignette has a segment-shuffling section. NEWS.md has an entry.
- [ ] AC7: `pkgdown::check_pkgdown()` passes. `devtools::document()` leaves no uncommitted diff in `man/` or `NAMESPACE`. `devtools::test()` passes. `devtools::check()` gives 0 errors, 0 warnings, and no NOTE that `main` does not also give.

## Coverage

- AC1 → T2
- AC2 → T2
- AC3 → T3
- AC4 → T2
- AC5 → T4
- AC6 → T1, T5
- AC7 → T6

## Tasks

- [x] T1: Write `cairn/references/tschacher2020.md` from the source-note template. The source is the SUSY 0.1.0 code on the CRAN GitHub mirror (`R/susy.R`, mirror commit `02590aea0877f09411584827cecaa623e46b6e4f`). A copy is at `sources/susy-0.1.0-susy.R`. Extract the segment count (`round(size/range - 0.499999)`), the dropped tail, the two segment-size rules, and the all-pairs surrogate set. State where bsync departs, including the seconds-to-rows check and the segment count on complete cases. Add the INDEX.md line.
- [x] T2: Write the AC1, AC2, and AC4 tests first. Put the plain reference in the test file with an oracle provenance header (DESIGN §13). Then implement `generate_surrogate_segment()` next to the other within-dyad generators in `R/surrogate_generation.R`. For small k, enumerate every order with no fixed segment and sample from that list. For large k, reject repeated draws.
- [ ] T3: Write the AC3 calibration test. Pick the series length, segment size, and lag so that the test runs in under a minute. Record the observed rejection rate in the work log.
- [ ] T4: Add `"segment"` to the `surrogate_method` choices in `synchrony_multiverse()` and to the `autotune_wcc()` docs. Generate one matrix per distinct window size and look it up per cell. Skip cells whose window is above half the series length, with the existing skip path (`R/multiverse.R`, about line 303). Write the AC5 tests.
- [ ] T5: Write roxygen with `@md`, `@references`, and `@seealso` links both ways with the other generators. Add the `_pkgdown.yml` row, the vignette section, the NEWS entry, and the WORDLIST words. Update the method count in section 1 of the vignette. Update DESIGN §2, §6, §13 (oracle records), and §14 #8. Run `devtools::document()`.
- [ ] T6: Run `devtools::test()`, `spelling::spell_check_package()`, `pkgdown::check_pkgdown()`, and `devtools::check()`. Run `air format` only on files that were air-clean before the edit.

## Work log

- 2026-10-06: created by /milestone-plan. It takes the segment-shuffling half of the "Expanded surrogate generators" candidate row. M014 takes the IAAFT half. It depends on M014 because both edit the same `surrogate_method` choices and vignette.
- 2026-10-06: question set: can the agent fetch the primary sources itself — yes. The agent read the SUSY 0.1.0 source through the GitHub API and copied it to the shelf.
- 2026-10-06: question set: offer the new generators in `synchrony_multiverse()` and `autotune_wcc()` — yes, both. The segment size comes from each cell's window size.
- 2026-10-06: plan gate chose the reordered series over SUSY-exact pairing, because it works with all four wrappers and the multiverse; falsified by a user need to reproduce SUSY's numbers exactly.
- 2026-10-06: plan chose to keep the tail in place over dropping it as SUSY does, because each wrapper needs `length(y)` rows; falsified by a calibration run where tail windows inflate the rejection rate.
- 2026-10-06: plan chose distinct orders without replacement over independent draws, because draws with replacement never include the original order and reject too often at small k (phipson2010, pp. 6-7); falsified by a calibration run at k = 3 to 5 with a rejection rate above 0.09.
- 2026-10-06: criteria audit (full mode, fresh Opus reader): 6 findings, all fixed. Draws with replacement gave invalid p-values at small k, so orders are now distinct (Scope, AC3, AC4). AC3 now fixes the AR coefficient, k, and lag_max. AC5 skips a too-long window instead of aborting, and names the lag and increment. Scope names SUSY's seconds-to-rows check. AC1 adds a hand-derived SUSY segment cut. AC7 no longer binds the Review record.
- 2026-10-06: oracle plan: closed-form (plain reimplementation, AC2), invariant (AC1, AC5), and simulation-coverage (AC3). `test-external-oracle.R` already pins the SUSY kernel, and this generator only reorders `y`.
- 2026-10-07: implement started on branch m015-segment-shuffle-surrogates, cut from origin/main at 9f82936.
- 2026-10-07: T1 done. `cairn/references/tschacher2020.md` written from SUSY 0.1.0 `R/susy.R`, `man/susy.Rd`, and `DESCRIPTION` (the DOI is read from `DESCRIPTION`). INDEX line added.
- 2026-10-07: T2 choice: up to 8 segments (14833 orders at most), the generator lists every order with no fixed segment in lexicographic order. One `sample.int()` call then draws distinct rows. Above 8, it draws `sample.int(k)` and rejects a fixed segment or a repeat. Above 8 segments, `n_surrogates` at or above the order count aborts, because listing 9! or more orders is too large.
- 2026-10-07: T2 choice: the tail and all-orders messages are `cli_inform()` with the classes `bsync_segment_tail` and `bsync_segment_all_orders`. If the tail is empty, no tail message fires.
- 2026-10-07: T2 done. `generate_surrogate_segment()` and `tests/testthat/test-surrogate-segment.R` (AC1, AC2, AC4) added, with the `_pkgdown.yml` row. Planted defects (a fixed segment allowed, the tail moved) turned 5 and 36 expectations red. Full suite: 0 failures.
- 2026-10-07: T3 finding: the planned design is liberal. Size runs used 1000 independent AR(1) pairs (phi 0.7), 19 surrogates, `wcc_surrogate()` with window 32 and lag 4, at p <= .05. With orders that keep no segment in place and segment = window = increment, the rate was 0.090 at 8 segments, 0.118 at 5, and 0.073 at 20. With a 3-sample increment (the multiverse default of 10%), it was 0.104. Phase surrogates gave 0.055 to 0.059.
- 2026-10-07: T3 finding: two causes. (1) Orders that keep no segment in place are not a group with the identity, so the surrogates resemble each other more than the observed series. Drawing any order except the identity gave 0.049 at 5 segments and 0.068 at 8. (2) Lagged windows cross segment boundaries, so the observed series keeps continuity that the surrogates lose. With segment = window + 2 * lag_max = increment, every lagged window stays in one segment. Both fixes together gave 0.049 at 8 segments and 0.043 at 5, and 0.060 with a 3-sample increment.
- 2026-10-07: stop: the fix changes the Goal (no segment in place) and the multiverse segment size the question set chose (the window size). The user decides between a re-plan, an escalation, and the plan as written.
- 2026-10-07: user chose escalation to a Fable review through /milestone-brief before any re-plan.
- 2026-10-07: blocked on RB01 (`cairn/reviews/RB01-segment-surrogate-size.md`). The brief commit goes on the milestone branch, not main, because main still carries the M015 file at `planned`. The calibration script is inline in the brief, because the session scratchpad is not in the repo.

## Decisions

## Review

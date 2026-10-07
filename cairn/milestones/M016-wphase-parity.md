# M016: wphase parity: selectable statistic and workflow article

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Resolves:** —
- **Surface tier:** user-facing — adds an exported argument and a pkgdown article
- **Branch/PR:** m016-wphase-parity

## Goal

Bring `wphase()` to the level of the other estimators: a selectable aggregate statistic that its surrogate test matches, and a workflow article.

## Scope

**In:**
- A `statistic = c("mean_plv", "peak")` argument on `wphase()` and `wphase_surrogate()`, default `"mean_plv"` (current behavior). `"peak"` takes the largest PLV across lags in each window, then averages those per-window values. This is the definition `wcc_aggregate()` uses for wcc's `"peak"` (`R/wcc.R:386`), applied to PLV. The surrogate aggregate uses the same helper as the observed one (CLAUDE.md Invariant 2).
- The chosen statistic in `settings$statistic`, in `names(aggregate)`, in `glance()`, and in the `print()` label of `wphase_res` and `wphase_surr`.
- A new article `vignettes/wphase-workflow.Rmd`, built like the wcc/wdtw/wgranger workflow articles, with band-pass guidance before phase extraction. The band-pass uses `gsignal::butter()` and `gsignal::filtfilt()` directly (gsignal is already in Imports).
- Roxygen for the new argument, NEWS.md, the DESIGN.md §9 defaults row, `_pkgdown.yml` articles list.

**Out:**
- wphase in `synchrony_multiverse()` / `autotune_wcc()`: not requested. It gets a new candidate row if a user asks.
- A band-pass method in `smooth_signal()`: the article uses gsignal directly (see work log).
- wphase segment surrogates on the phase series: stays the `[low]` "wphase segment surrogates" candidate row.
- The window-count fix for `window_increment > 1`: M017.

## Acceptance criteria

- [x] AC1: `wphase()` and `wphase_surrogate()` accept `statistic = "mean_plv"` (the default) and `statistic = "peak"`. Tests in `tests/testthat/test-wphase.R` compute the `"peak"` value with an explicit loop (per-window maximum PLV across lags, then the mean over windows) and find `wphase(..., statistic = "peak")$aggregate` equal to it within 1e-12 in two cases:
      (a) a hand-built phase pair with default increments and no `time`,
      (b) `window_increment = 2`, `lag_increment = 3`, and a `time` vector with one tied pair of timestamps.
      With no `statistic` argument, `wphase()$aggregate` is `identical()` to `base::mean(results_df$plv, na.rm = TRUE)` in both cases.
- [x] AC2: For each `statistic` value and each AC1 case, a test runs `wphase_surrogate()` on a 3-column surrogate matrix. It finds each `surrogate_z[j]` equal to `wphase(x, y_surrogates[, j], statistic = )$aggregate` within 1e-12, and `observed_z` for `"peak"` equal to the AC1 explicit-loop value.
- [x] AC3: A test finds that an unknown `statistic` value makes `wphase()` and `wphase_surrogate()` stop with an error that matches "should be one of".
- [x] AC4: For each `statistic` value, a test finds the value in `settings$statistic`, `names(aggregate)`, and `glance()$statistic` of a `wphase_res`. Asserted with `expect_message()` (LESSONS M008), `print()` of a `wphase_res` shows "Mean PLV" or "Mean Peak PLV", and `print()` of a `wphase_surr` shows "Observed Mean PLV" / "Average Null Mean PLV" or "Observed Mean Peak PLV" / "Average Null Mean Peak PLV".
- [x] AC5: `vignettes/wphase-workflow.Rmd` exists and is listed under articles in `_pkgdown.yml` next to `wgranger-workflow`. Its section headings cover these 8 topics:
      (1) what windowed phase synchrony measures and its narrowband assumption,
      (2) band-pass filtering before phase extraction,
      (3) `wphase()` on the filtered pair,
      (4) a surrogate test with `wphase_surrogate()`,
      (5) `statistic = "peak"`,
      (6) `pick_optima()`,
      (7) plots,
      (8) a See also list. Every result its prose states (a PLV, a p-value, a window count) comes from inline R code or chunk output.
- [x] AC6: NEWS.md has an entry for the `statistic` argument and the new article, with no milestone numbers.
- [x] AC7: `devtools::document()` produces no diff, `pkgdown::check_pkgdown()` passes, and `devtools::check()` gives 0 errors and 0 warnings (any NOTE justified in the Review section).

## Coverage

- AC1 → T1, T2
- AC2 → T1, T2
- AC3 → T1, T2
- AC4 → T3
- AC5 → T4
- AC6 → T5
- AC7 → T5

## Tasks

- [x] T1: Write the tests for AC1–AC3 in `tests/testthat/test-wphase.R` first. The hand-built pair uses phases with a known PLV per window and lag (for example a constant phase offset in some windows and scattered phases in others), so the explicit loop has a known answer.
- [x] T2: Add `statistic` to `wphase()` (`R/wphase.R:65`) and `wphase_surrogate()` (`R/surrogate_analysis.R`), with `match.arg()`. Extend `wphase_aggregate()` (`R/wphase.R:211`) to take the statistic and the window ids, as `wcc_aggregate()` does, and call it from both paths. Take the window ids from the grid positions before any `time` mapping, so that tied timestamps group windows the same way on the observed and surrogate sides. Keep the `"mean_plv"` path as `base::mean(plv, na.rm = TRUE)`. Document the argument with the reason for the default and the source of the `"peak"` definition.
- [x] T3: Carry the statistic into `settings`, `names(aggregate)`, and the `print()` labels of `print.wphase_res` and `print.wphase_surr`. Write the AC4 tests.
- [x] T4: Write `vignettes/wphase-workflow.Rmd` from the structure of `vignettes/wgranger-workflow.Rmd` and `wdtw-workflow.Rmd`. Show the effect of the band-pass on mean PLV with computed values, keep the surrogate chunk fast enough for `check()`, and do not set a `future` multisession plan (LESSONS M011). Add it to `_pkgdown.yml` and link it from the See also lists of the other workflow articles and from `wphase()` `@seealso`.
- [x] T5: Write the NEWS.md entry. Add the DESIGN.md §9 row "wphase aggregate statistic: selectable, default `"mean_plv"`". Run `spelling::spell_check_package()` and update `inst/WORDLIST` (LESSONS M011). Run `air format` only on files that were air-clean before the edit (LESSONS M010). Run `devtools::document()`, `devtools::test()`, `pkgdown::check_pkgdown()`, and `devtools::check()`.

## Work log

- 2026-10-07: created by /milestone-plan. Absorbs the candidate rows "wphase workflow vignette" and "Selectable wphase aggregate statistic" (both from M008 Out).
- 2026-10-07: question set: what to plan next — wphase parity.
- 2026-10-07: plan gate chose band-pass with `gsignal::butter()` + `gsignal::filtfilt()` in the article over a new band-pass method in `smooth_signal()`, because no one asked for new API and gsignal is already imported; falsified by users needing band-pass often enough that the article's code is copied into analyses.
- 2026-10-07: plan gate chose `"peak"` as the second statistic over other PLV summaries (median, lag-0 only) because it matches wcc's `statistic` values, so one argument means the same thing across estimators; falsified by a published PLV convention that defines a different per-window summary.
- 2026-10-07: criteria audit (full mode, fresh Opus reader) returned 6 findings on M016, all fixed: AC1 pins the default path to `mean(results_df$plv, na.rm = TRUE)` instead of "the value on main" (M017 moves fixtures with increments above 1); AC1/AC2 add a case with increments above 1 and a tied `time`, and T2 takes window ids from grid positions; AC2 drops the trivially true `observed_z` clause for an explicit-loop check; AC3 matches "should be one of" only; AC4 names the print labels; AC5 limits the computed-number rule to results. The tied-time issue in `wcc()` "peak" became a `[low]` candidate row.
- 2026-10-07: T1 done. Added 3 tests to `test-wphase.R` (peak oracle loop, matched null, unknown statistic). All 3 fail on main with "unused argument (statistic = ...)". Case (b) uses lag_max = 5 so the lags (-5, -2, 1, 4) are asymmetric, and n = 201 so its window count is the same before and after M017.
- 2026-10-07: T2 done. `create_wphase_df()` now returns `list(results_df, window_id)`, with `window_id` = grid positions before the time mapping. `wphase_aggregate(plv, window_id, statistic)` serves both paths. Added `@md` to the `wphase()` block (the per-block pattern of M009/M011/M013). Planted-defect checks: grouping by time-mapped `i` failed the AC1 and AC2 tests (1 and 4 failures), a surrogate side fixed to `"mean_plv"` failed AC2 (6 failures). The first tie (starts 7 and 9) hit no window start and missed the plant, so the tie moved to starts 6 and 8. Full suite: 313 tests, 0 failed.
- 2026-10-07: T3 done. One internal `wphase_agg_label()` gives the label to both print methods, and both use the `setNames()` card form of the wcc prints. `glance()` already copies `settings$statistic` (`R/tidy.R:94`), so it needed no change. The AC4 test also checks that "Mean PLV" is absent from the "peak" card.
- 2026-10-07: T4 done. `vignettes/wphase-workflow.Rmd` simulates a 0.5 Hz rhythm with shared tempo wander (LESSONS M008), independent AR(0.98) drift, and a lead shifting from +0.5 s to -0.5 s, then band-passes 0.3-0.7 Hz. It renders in about 4 s with 1000 circular surrogates (`window_increment = 5`). Implementation choice: the article has no `leadership_asymmetry()` section, unlike the wdtw/wcc articles. A prototype showed the PLV optima correlate 0.50 with the simulated lead (relative phase at lag 0: 0.98), and the LAI pointed the wrong way in the third quarter. The article states both correlations from code and advises reading lead-lag from the relative phase. The wcc article had no See also section, so one was added for the link.
- 2026-10-07: T5 done. NEWS entry under "Windowed phase synchrony"; DESIGN §9 row added. The spell check flagged the Greek letters in the article's lead formula, so that sentence now uses words and `2 * pi * 0.5` (no WORDLIST change). `document()` no diff, `check_pkgdown()` no problems, `devtools::check()` 0 errors / 0 warnings / 0 notes. `air format` not run: `R/wphase.R`, `R/surrogate_analysis.R`, and `test-wphase.R` were not air-clean on main (LESSONS M010).
- 2026-10-07: claim audit: 61 claims read, 5 corrected — R/surrogate_analysis.R, vignettes/wphase-workflow.Rmd, tests/testthat/test-wphase.R (statistic param and Section 6 now say the test is of the reported aggregate; "7.5 cycles"; circular surrogates take `lag_max = 10` so no shift lies within the lag range, p now 0.001 for both statistics; case (b) comment states the grid arithmetic). The same reader re-read all 5 and found them correct. `devtools::check()` after the fixes: 0 errors / 0 warnings / 0 notes.
- 2026-10-07: all tasks done; status set to review.

## Decisions

## Review

Evidence gathered 2026-10-07 at 2e7faad, branch current with origin/main (8abadee).

- AC1: `test_file("test-wphase.R")` gives 16 tests, 0 failed, 0 errors, 0 skipped. Test at line 361 runs both cases from `stat_cases()`: (a) the hand-built pair (default increments, no `time`), and (b) `window_increment = 2`, `lag_increment = 3`, and `time` with starts 6 and 8 tied. It checks `"peak"` against `peak_plv_loop()` with `expect_lt(abs(diff), 1e-12)` and the default with `expect_identical(..., base::mean(results_df$plv, na.rm = TRUE))`. In T2, grouping by the time-mapped `i` failed this test.
- AC2: the same run passes the test at line 375. For both cases, both statistics, and a 3-column `generate_surrogate_circular()` matrix, it checks each `surrogate_z[j]` against `wphase(x, y_surrogates[, j], statistic = )$aggregate` within 1e-12, and `observed_z` for `"peak"` against `peak_plv_loop()`. In T2, a surrogate side fixed to `"mean_plv"` gave 6 failures.
- AC3: the same run passes the test at line 423. `wphase(..., statistic = "median")` and `wphase_surrogate(..., statistic = "median")` each raise an error matching "should be one of". Before T2, the same test failed with "unused argument (statistic = ...)".
- AC4: the same run passes the test at line 398. For `"mean_plv"` and `"peak"` it checks `settings$statistic`, `names(aggregate)`, and `glance()$statistic` with `expect_identical`. With `expect_message(..., fixed = TRUE)` it finds "Mean PLV" / "Mean Peak PLV" in the `wphase_res` print and "Observed …" / "Average Null …" with the same labels in the `wphase_surr` print. It also finds no "Mean PLV" in the `"peak"` card.
- AC5: `vignettes/wphase-workflow.Rmd` exists, and `_pkgdown.yml:32` lists `wphase-workflow` directly after `wgranger-workflow`. Headings (`grep '^##'`) map to the 8 topics as follows. (1) §1 "What is Windowed Phase Synchrony?" with §1.1 "The Narrowband Assumption". (2) §3 "Band-Pass Filtering Before Phase Extraction". (3) §4 "Calculating Windowed Phase Synchrony on the Filtered Pair". (4) §5 "Surrogate Testing for Significance", which calls `wphase_surrogate()`. (5) §6 "The Peak Statistic". (6) §7 "Optima Extraction". (7) §8 "Visualizing the Results". (8) "See also". An awk sweep of prose lines (outside chunks) that carry decimals or "p =" found every result as inline R. The remaining literals are simulation or setup parameters (0.5 Hz, 7.5 cycles, the simulated ±0.5 s lead). The article built inside `devtools::check()`.
- AC6: NEWS.md "Windowed phase synchrony" (lines 3-15) has one bullet for the `statistic` argument and one for `vignette("wphase-workflow")`. `grep -n "M0" ` over those lines finds no milestone number, and `test-doc-hygiene.R` passed inside `check()`.
- AC7: `devtools::document()` left `git status` clean. `pkgdown::check_pkgdown()` reported "No problems found". `devtools::check()` at 2e7faad gave 0 errors, 0 warnings, and 0 notes.

Consistency gate: `cairn_validate.py` exit 0 ("all checks passed"). No principle changed, so `cairn_impact` was skipped. Profile slot results: `document()` gives no diff. Generated files are not hand-edited (man/ comes from `document()`). The branch does not touch README.Rmd or README.md. `check_pkgdown()` passes. The NEWS entry is present. There are no new top-level files. `check()` gives 0/0/0.

Independent review. spawned: diff-bug, blame-history, prior-review. The 22 findings were settled as 12 fixed, 2 sent to candidate rows, and 8 rejected. After the fixes, `devtools::test()` gave 314 tests and 0 failed, `devtools::check()` gave 0/0/0, and the spell check was clean.

- prior-review #1: the new article is missing from the article lists in `vignettes/bsync.Rmd`, README, and `choosing-parameters.Rmd` — fix now: added to the bsync.Rmd and README "Where to go next" tables and to the choosing-parameters See also list, and README.md rebuilt (one added line), fixed d74beab. The getting-started estimator table had no wphase column before this branch, so that part became the candidate row "Getting-started estimator table lacks wphase".
- prior-review #2: `@md` on the `wphase()` block re-renders its prose as markdown — reject (false as a defect): the diff-bug lens read `man/wphase.Rd` and found that the change fixes literal backticks and `[fn()]` that rendered raw on main. This follows the per-block pattern of M009/M011/M013.
- prior-review #3: qualitative claims in the article (peak larger than mean for the surrogates, heatmap "high across most lags") are not computed — fix now: the peak-versus-mean comparison is now an inline count (1000 of 1000 surrogates), and the heatmap sentence states only what the inline medians show, fixed d74beab.
- diff-bug #1: with tied `time`, `print()` / `glance()` / `pick_optima()` count windows by the time-mapped `i` (79) while the "peak" aggregate uses grid positions (80) — follow-up: the undercount was there before this branch. Absorbed into the candidate row "Tied timestamps merge windows".
- diff-bug #2: the article's "drift decides the PLV, not the coupling" compares unlike pairs — fix now: it now says that two series that share nothing score higher than a coupled pair, so unfiltered PLV does not measure coupling. Fixed d74beab.
- diff-bug #3: the article does not mention optima at the lag boundary — fix now: an inline count of optima at ±10 samples (44 of 206) is now stated without a causal claim, fixed d74beab.
- diff-bug #4: the NEWS claim that relative phase tracks lead–lag better turns one simulation into a rule — fix now: it now says "In its simulated example", fixed d74beab.
- diff-bug #5: "rMEA best-lag convention (Boker et al., 2002)" mixes sources — fix now for `wphase()`: the docs now call it the best-lag definition that `wcc()` uses, with no citation, fixed d74beab. The same wording in `wcc()` existed before this branch and became the candidate row "wcc "peak" citation".
- diff-bug #6: "high across most lags in every window" is not computed — fix now: same fix as prior-review #3, fixed d74beab.
- diff-bug #7: "second-order Butterworth" — fix now: it now says `butter()` with `n = 2` gives a 4th-order band-pass, fixed d74beab.
- diff-bug #8: `wphase_aggregate()` treats any non-`"mean_plv"` value as `"peak"` — reject (planned change): T2 specified the shape of `wcc_aggregate()`, and both exported entry points validate with `match.arg()`.
- diff-bug #9: the AC4 test leaks the extra cli lines of the surrogate print into the test output — reject (style): output noise only, and AC4 requires `expect_message()`.
- diff-bug #10: the AC1 default check repeats the implementation's expression — reject (planned change): AC1 specifies exactly this `identical()` check, and the "peak" path has the explicit-loop oracle.
- diff-bug #11: the history wording "the value `wphase()` reported before this argument existed" in the roxygen — fix now: removed, fixed d74beab.
- blame-history #1: the article cautions against `leadership_asymmetry()` with wphase optima, which M008 made supported — fix now: the article now says the function accepts wphase optima and advises checking them against the relative phase first, fixed d74beab.
- blame-history #2: "as `wcc(statistic = "peak")` does" is not exact under tied `time` — follow-up: the difference is wcc's tied-time grouping. Absorbed into the candidate row "Tied timestamps merge windows".
- blame-history #3: the peak oracle never exercises the window that M017 will add — reject (false): `peak_plv_loop()` enumerates every fitting start on its own (`i + ws - 1 + lag_max <= n`), so it stays correct after M017. n = 201 only keeps the current grid equal to it.
- blame-history #4: temporal wording in the roxygen — fix now: same fix as diff-bug #11, fixed d74beab.
- blame-history #5: `wphase_aggregate()` falls into "peak" for any other value — reject (planned change): same as diff-bug #8.
- blame-history #6: `@md` re-renders the whole Rd — reject (false as a defect): same as prior-review #2.
- blame-history #7: NEWS overgeneralizes from one simulation — fix now: same fix as diff-bug #4, fixed d74beab.
- blame-history #8: the shared `wphase_agg_label()` helper differs from wcc's inline labels — reject (style).

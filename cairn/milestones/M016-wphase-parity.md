# M016: wphase parity: selectable statistic and workflow article

- **Status:** in-progress
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

- [ ] AC1: `wphase()` and `wphase_surrogate()` accept `statistic = "mean_plv"` (the default) and `statistic = "peak"`. Tests in `tests/testthat/test-wphase.R` compute the `"peak"` value with an explicit loop (per-window maximum PLV across lags, then the mean over windows) and find `wphase(..., statistic = "peak")$aggregate` equal to it within 1e-12 in two cases:
      (a) a hand-built phase pair with default increments and no `time`,
      (b) `window_increment = 2`, `lag_increment = 3`, and a `time` vector with one tied pair of timestamps.
      With no `statistic` argument, `wphase()$aggregate` is `identical()` to `base::mean(results_df$plv, na.rm = TRUE)` in both cases.
- [ ] AC2: For each `statistic` value and each AC1 case, a test runs `wphase_surrogate()` on a 3-column surrogate matrix. It finds each `surrogate_z[j]` equal to `wphase(x, y_surrogates[, j], statistic = )$aggregate` within 1e-12, and `observed_z` for `"peak"` equal to the AC1 explicit-loop value.
- [ ] AC3: A test finds that an unknown `statistic` value makes `wphase()` and `wphase_surrogate()` stop with an error that matches "should be one of".
- [ ] AC4: For each `statistic` value, a test finds the value in `settings$statistic`, `names(aggregate)`, and `glance()$statistic` of a `wphase_res`. Asserted with `expect_message()` (LESSONS M008), `print()` of a `wphase_res` shows "Mean PLV" or "Mean Peak PLV", and `print()` of a `wphase_surr` shows "Observed Mean PLV" / "Average Null Mean PLV" or "Observed Mean Peak PLV" / "Average Null Mean Peak PLV".
- [ ] AC5: `vignettes/wphase-workflow.Rmd` exists and is listed under articles in `_pkgdown.yml` next to `wgranger-workflow`. Its section headings cover these 8 topics:
      (1) what windowed phase synchrony measures and its narrowband assumption,
      (2) band-pass filtering before phase extraction,
      (3) `wphase()` on the filtered pair,
      (4) a surrogate test with `wphase_surrogate()`,
      (5) `statistic = "peak"`,
      (6) `pick_optima()`,
      (7) plots,
      (8) a See also list. Every result its prose states (a PLV, a p-value, a window count) comes from inline R code or chunk output.
- [ ] AC6: NEWS.md has an entry for the `statistic` argument and the new article, with no milestone numbers.
- [ ] AC7: `devtools::document()` produces no diff, `pkgdown::check_pkgdown()` passes, and `devtools::check()` gives 0 errors and 0 warnings (any NOTE justified in the Review section).

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
- [ ] T5: Write the NEWS.md entry. Add the DESIGN.md §9 row "wphase aggregate statistic: selectable, default `"mean_plv"`". Run `spelling::spell_check_package()` and update `inst/WORDLIST` (LESSONS M011). Run `air format` only on files that were air-clean before the edit (LESSONS M010). Run `devtools::document()`, `devtools::test()`, `pkgdown::check_pkgdown()`, and `devtools::check()`.

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

## Decisions

## Review

# M013: Granger multiverse summary and plot for both directions

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Resolves:** —
- **Surface tier:** user-facing — changes exported print, summary, glance, and plot methods
- **Branch/PR:** m013-granger-multiverse-direction

## Goal

Let users read the y -> x direction of a Granger `synchrony_multiverse()` result through the print, summary, glance, and plot methods that now report only x -> y.

## Scope

**In:**
- `synchrony_multiverse(estimator = "wgranger")` also stores `$robustness_yx`. It uses the grid's `es_yx` and `p_yx` columns and has the same fields as `$robustness`. `$robustness` stays x -> y.
- A `direction = c("xy", "yx")` argument, default `"xy"`, on `print()`, `summary()`, `glance()`, and `plot()` for `bsync_multiverse`. The value `"yx"` on a WCC or WDTW result is an error.
- For Granger results, `glance()` gets a `direction` column. `print()`, `summary()`, and the plot title name the direction that they report. WCC and WDTW output does not change.
- A test for the "No valid cells to plot" error in both directions. This absorbs the candidate row.
- The `choosing-parameters` vignette, the `R/multiverse.R` header comment, and one test title name an `es_xy` grid column that does not exist. The x -> y columns are `es` and `p`. Fix them, and show the y -> x call in the vignette.
- A NEWS entry.

**Out:**
- A view that shows both directions at once (two panels or two curves). The plan gate chose the `direction` argument instead.
- `autotune_wcc()` and `select_specification()`. They are WCC only, so they have no Granger direction.
- Other multiverse and wphase candidates stay as ROADMAP candidate rows.

## Acceptance criteria

- [x] AC1: A Granger `synchrony_multiverse()` result has `$robustness_yx` with the same seven field names as `$robustness`. The test fixture has `n_significant` and `median_es` values that differ between the two directions. The test computes all seven fields directly from the grid's `es_yx` and `p_yx` columns, with significance as `p <= .05`, and shows that `$robustness_yx` holds these values. The same test shows that `$robustness` holds the same computation over `es` and `p`. Tests show that WCC and WDTW results have no `$robustness_yx` element.
- [x] AC2: `glance(direction = "yx")` on a Granger result returns the `$robustness_yx` values in its summary columns and `"yx"` in a `direction` column. With the default, it returns the `$robustness` values and `"xy"`. For WCC and WDTW results, a test shows that the `glance()` column names are exactly `estimator, n_cells, n_valid, n_significant, pct_significant, median_es, iqr_es, sign_consistent, n_surrogates`.
- [x] AC3: On a Granger result, `print()` and `summary()` name the direction they report, for `"xy"` and for `"yx"`. With `"yx"`, they report the `$robustness_yx` counts, and `summary()` reports the `es_yx` range and counts skipped or NA cells from `es_yx`. Tests assert each of these with `expect_message()` matchers on the direction line, the significant-count line, and the ES-range line, on the AC1 fixture.
- [x] AC4: `plot(direction = "yx")` on a Granger result ranks cells by `es_yx` and marks significance from `p_yx`. A test on the plot-data helper shows this on a fixture where the two directions rank cells differently and at least one cell is significant in one direction only. A test asserts the plot title string for `"xy"` and `"yx"` on a Granger result, and for a WCC result. A y -> x vdiffr snapshot is committed after a visual check. `git diff main -- tests/testthat/_snaps/multiverse/multiverse-spec-curve.svg` is empty.
- [x] AC5: For each of `print()`, `summary()`, `glance()`, and `plot()`, `direction = "yx"` on a WCC result gives an error that names the estimator. A test shows the same error for a WDTW result on at least one method. For each of the four methods, an unknown `direction` value gives an error. Each test asserts the message with a regexp.
- [x] AC6: Tests reach the "No valid cells to plot" error on two hand-built grids: one with `es` non-NA and `es_yx` all NA (with `direction = "yx"`), and one the reverse (with `direction = "xy"`). Each test asserts the message with a regexp.
- [x] AC7: Roxygen documents `direction` on the four methods and `$robustness_yx` in the return value of `synchrony_multiverse()`. `grep -rn "es_xy" vignettes/ tests/` returns no lines. The `choosing-parameters` vignette states that `es` and `p` are x -> y and calls `glance(mv_g, direction = "yx")`. Lines 1-15 of `R/multiverse.R` name `es` and `p` as the x -> y columns. NEWS.md has an entry for this change. `devtools::document()` leaves no uncommitted diff in `man/` or `NAMESPACE`. `devtools::test()` passes. `devtools::check()` gives 0 errors and 0 warnings, and each NOTE is triaged in the Review section.

## Coverage

- AC1 → T1
- AC2 → T2
- AC3 → T3
- AC4 → T4
- AC5 → T2, T3, T4
- AC6 → T4
- AC7 → T5, T6

## Tasks

- [x] T1: In `synchrony_multiverse()` (`R/multiverse.R:426`), compute `multiverse_robustness(n_cells, es_yx_vec, p_yx_vec, skipped_vec)` for Granger and store it as `$robustness_yx`. Add one internal helper that matches `direction` with `rlang::arg_match()`, stops for `"yx"` on a non-Granger result, and returns the matching robustness list and grid columns. Write the AC1 tests first. Pick a seed and series where `n_significant` and `median_es` differ between directions, and where at least one cell is significant in one direction only.
- [x] T2: Add `direction` to `glance.bsync_multiverse()` (`R/tidy.R:154`). Add the `direction` column for Granger only. Write the tests for AC2 and for the `glance()` part of AC5.
- [x] T3: Add `direction` to `print.bsync_multiverse()` (`R/multiverse.R:568`) and `summary.bsync_multiverse()` (`R/multiverse.R:590`). Add a direction line for Granger. For `"yx"`, `summary()` reads `es_yx` for the range and the skipped count. Write the tests for AC3 and for their parts of AC5.
- [x] T4: Add `direction` to `plot.bsync_multiverse()` and `multiverse_plot_data()` (`R/multiverse_plot.R:12`, `R/multiverse_plot.R:42`). Move the title text (`R/multiverse_plot.R:73`) into an internal helper. It reads the chosen robustness list and names the direction for Granger only, so the WCC title stays the same. Write the tests for AC4, AC6, and the plot part of AC5. Record the new vdiffr snapshot and look at it before the commit.
- [x] T5: Write roxygen for `direction` and `$robustness_yx`, with `@md` per the roxygen-markdown lesson. Fix the comment at `R/multiverse.R:13` and the test title at `tests/testthat/test-multiverse.R:293`. Fix `vignettes/choosing-parameters.Rmd:151-166` and add the y -> x `glance()` call. Add the NEWS.md entry. Run `devtools::document()`.
- [x] T6: Run `devtools::test()`, `spelling::spell_check_package()`, and `devtools::check()`. Run `air format` only on files that were air-clean before the edit.

## Work log

- 2026-10-06: created by /milestone-plan. The request named no work. The question set offered four candidate rows, and the user chose the Granger multiverse direction gap.
- 2026-10-06: collision check: absorbs the candidate rows "Multiverse summary and plot ignore the y -> x Granger direction" and "Test the all-NA error in `plot.bsync_multiverse()`". No D-entry, planned milestone, or archive overlaps. Inbox sweep found 0 open issues and 0 open PRs.
- 2026-10-06: plan gate chose a `direction` argument (default `"xy"`) over always showing both directions, because the default output stays the same. Falsified by users who need both directions in one call, such as one `glance()` row per direction.
- 2026-10-06: plan decided that `glance()` gets a `direction` column for Granger only, so WCC and WDTW columns do not change. The `es_xy` doc fix is in scope because the same vignette section shows the new call.
- 2026-10-06: criteria audit (full mode, fresh Opus reader) returned 13 findings, and all 13 are fixed. Fixtures must differ in `n_significant`, `median_es`, ranking, and one significance call. All seven robustness fields are checked, and WDTW cases are added. The glance column list is written out. The plot title moves to a tested helper, and snapshot checks use `git diff`. AC6 uses hand-built grids. The `R/multiverse.R` header and a test title are named, and NOTEs are triaged. Proportionality: no finding.
- 2026-10-06: implement started on branch `m013-granger-multiverse-direction`.
- 2026-10-06: T1 done. `$robustness_yx` is set for Granger, and the `multiverse_direction()` helper matches `direction` and stops for `"yx"` on WCC or WDTW. The fixture is an AR(1) x driving y at lag 2, with no y -> x path: 8 cells, 6 significant x -> y and 0 y -> x. The header comment fix (part of T5) landed here. The touched R files were not air-clean on main, so only the new test file is air-formatted. Suite: 1921 pass, 0 fail.
- 2026-10-06: T2 done. `glance()` takes `direction` and adds a `direction` column for Granger only, with `tibble::add_column()`. The block gets `@md`. Suite: 1931 pass, 0 fail.
- 2026-10-06: T3 done. `print()` and `summary()` take `direction`, and print a "Direction" line for Granger only. `summary()` reads the chosen ES column for the range and the skipped/NA count. With the T3 code stashed, 10 of the new expectations failed. Suite: 1949 pass, 0 fail.
- 2026-10-06: T4 done. `plot()` takes `direction` as its last named argument, after `top_frac`, so positional calls keep working. `multiverse_plot_data()` copies the chosen direction into `es`/`p`, and `multiverse_plot_title()` adds "(x -> y)" or "(y -> x)" for Granger only. The y -> x snapshot `multiverse-granger-yx.svg` was rendered to PNG and checked by eye: 0 of 8 significant, all ES negative. `git diff main` on the old snapshot is empty. Suite: 1969 pass, 0 fail.
- 2026-10-06: T5 done. Roxygen documents `direction` on the four methods, which now carry `@md`, and documents `$robustness_yx` and all five y -> x grid columns in `synchrony_multiverse()`. That block stays without `@md`, as on main. The vignette and the test title are fixed, and the NEWS entry is added. `devtools::document()` rewrote 5 Rd files. `grep -rn "es_xy" vignettes/ tests/` returns no lines.
- 2026-10-06: T6 spell check flagged "summarise" 3 times in the new roxygen. Changed to American spelling, and the spell check is clean.
- 2026-10-06: claim audit: 60 claims read, 4 corrected — vignettes/choosing-parameters.Rmd, tests/testthat/test-multiverse.R, NEWS.md, R/multiverse.R.
- 2026-10-06: T6 done. `devtools::check()` at 78434ef gave 0 errors, 0 warnings, 0 notes. After the claim-audit commit, `devtools::test()` gives 1969 pass, 0 fail, and `devtools::document()` gives no diff. Status set to review.
- 2026-10-06: step-7 approval: m013-granger-multiverse-direction approved for merge

## Decisions

## Review

Evidence run 2026-10-06 on `m013-granger-multiverse-direction` at e6de332, which contains origin/main 57f510a. `devtools::test(filter = "multiverse")`: every test in `test-multiverse-direction.R` (15 tests) and `test-multiverse.R` passes, 0 failed, 0 skipped.

- AC1 evidence: "the fixture separates the two directions" (4 pass: n_significant 6 vs 0, median ES differs, order differs, a cell significant one way only). "a Granger result carries $robustness_yx built from es_yx and p_yx" (4 pass: names equal the 7 fields, all 7 equal `expected_robustness()` computed in the test with `p <= 0.05`, same for `$robustness` over `es`/`p`). "WCC and WDTW results have no $robustness_yx" (2 pass).
- AC2 evidence: "glance() reports the chosen Granger direction" (5 pass: default gives `direction = "xy"` and `$robustness` values, `"yx"` gives `"yx"` and `$robustness_yx` values, n_significant differs). "glance() columns of WCC and WDTW results do not change" (2 pass: `expect_named()` with the 9 listed columns for each).
- AC3 evidence: "print() names the Granger direction and reports its counts" (5 pass: `expect_message()` on "Direction: x -> y", "Direction: y -> x", and each direction's "Significant (p <= .05): b of n" line). "summary() reports the es_yx range and NA count for y -> x" (8 pass: direction lines, significant-count line, ES-range lines for `es` and `es_yx` that differ, "(including 3 skipped/NA)" for y -> x against "(including 0 skipped/NA)" for x -> y). All on the AC1 fixture.
- AC4 evidence: "plot data rank cells by es_yx and mark significance from p_yx" (6 pass, on the AC1 fixture, whose rankings differ and which has cells significant one way only). "the plot title names the Granger direction only" (7 pass: exact title strings for x -> y, y -> x, and WCC with no direction). "plot() draws the y -> x specification curve" (1 pass): `multiverse-granger-yx.svg` was rendered to PNG and checked by eye in implement (T4 work-log line). `git diff main -- tests/testthat/_snaps/multiverse/multiverse-spec-curve.svg | wc -l` gives 0.
- AC5 evidence: `expect_yx_refused()` asserts the regexp `needs a Granger result[\s\S]*estimator is "<estimator>"` with class `rlang_error`. It passes for `glance()` on WCC and on WDTW, and for `print()`, `summary()`, and `plot()` on WCC. `expect_bad_direction()` asserts "`direction` must be one of "xy" or "yx", not "zz"" and passes for all four methods.
- AC6 evidence: "plot() stops when the chosen direction has no valid cells" (4 pass). One grid has `es_yx` all NA: `plot(direction = "yx")` errors with "No valid cells to plot \(all ES values are NA\)", and the x -> y plot data still build. The other grid has `es` all NA: `plot(direction = "xy")` gives the same error, and the y -> x plot data still build.
- AC7 evidence: `\item{direction}` is in the print, summary, glance, and plot Rd files, and `man/synchrony_multiverse.Rd` names `robustness_yx` twice. `grep -rn "es_xy" vignettes/ tests/` exits 1 with no lines. `choosing-parameters.Rmd:153` says `es` and `p` are x -> y, and line 170 calls `glance(mv_g, direction = "yx")`. `R/multiverse.R:12-15` names observed/null_mean/null_sd/es/p as x -> y. NEWS.md line 3 opens the section "Both Granger directions in multiverse summaries and plots". `devtools::document()` writes nothing and leaves a clean tree. `devtools::test()` at e6de332 gives 1969 pass, 0 fail. `devtools::check()` at e6de332 gives 0 errors, 0 warnings, 0 notes, so no NOTE needs triage.

Consistency gate: `cairn_validate.py` passes all checks. No DESIGN principle changed, so `cairn_impact` is skipped. `devtools::document()` gives no diff. README.Rmd and README.md are not changed. `pkgdown::check_pkgdown()` finds no problems. NEWS.md has an entry. The branch adds no top-level files. `devtools::check()` gives 0/0/0.

spawned: diff-bug, blame-history, prior-review

- diff-bug #1: a Granger object without `$robustness_yx` (saved before this change, or built by hand) gives a `glance(direction = "yx")` with the 7 summary columns silently missing — fix now: `multiverse_direction()` recomputes the summary from `es_yx`/`p_yx` when it is missing, fixed e16657d.
- diff-bug #2: the same object gives unclear cli or math errors from `print()` and `plot()` — fix now, same fix as #1, fixed e16657d.
- diff-bug #3: a Granger object without `es_yx`/`p_yx` columns gives tibble warnings and a vctrs error, not a clear message — fix now: `multiverse_direction()` stops with a message when the columns are missing, fixed e16657d.
- diff-bug #4: `plot(res, "yx")` binds `"yx"` to `sig_color` — fix now: the `@param direction` text for `plot()` says to pass it by name. The argument stays last so positional color arguments keep working, fixed e16657d.
- diff-bug #5: with `es_yx` all NA, `summary(direction = "yx")` prints "ES range: [Inf, -Inf]" with a base R warning — fix now: the range line says "none" when no ES is computable, fixed e16657d (NEWS line in 7c10bb9).
- diff-bug #6: the `synchrony_multiverse()` roxygen block has no `@md`, so the new backtick text renders raw — follow-up: absorbed into the candidate row "Enable roxygen markdown package-wide". The block was without `@md` on main, and turning it on renders the whole block differently.
- diff-bug #7: the `$robustness` docs say cells not counted in `n_valid` "were skipped as too short", which is not the only reason for y -> x — fix now: the docs say cells with a computable ES in that direction, fixed e16657d.
- diff-bug #8: no y -> x cell in the fixture is significant, so `sign_consistent` for y -> x is only checked as NA — fix now: a test on a hand-built grid with significant `p_yx` cells checks all 7 fields, through the #1 recompute path, fixed e16657d.
- diff-bug #9: the plot title is tested through the helper, and only the snapshot shows that `plot()` passes it — reject, planned change: AC4 and T4 chose a tested title helper, and `multiverse_direction()` output is tested beside it.
- diff-bug #10: an object with no `settings$estimator` gets "estimator is ," in the refusal — fix now: the message falls back to "unknown", fixed e16657d.
- diff-bug #11: default Granger `glance()` and `print()` output gains a column and a line — reject, false: NEWS.md already says that Granger `print()` and `summary()` show a "Direction" line and that `glance()` has a `direction` column.
- blame-history #1: `$robustness_yx` reuses `multiverse_robustness()`, and the y -> x p-value source was not traced — noted: `R/surrogate_analysis.R:421` computes `p_val_yx` in the add-one form, so the summary follows the p-value rules on record.
- blame-history #2: says the candidate row "Test the all-NA error in `plot.bsync_multiverse()`" is still open — reject, false: the plan commit 57f510a removed it, and `grep -c "all-NA error" cairn/ROADMAP.md` gives 0.
- blame-history #3, #4, #6, #7, #8, #9, #10: no conflict found — noted.
- blame-history #5: Granger `glance()` gains a column, which could break strict column checks — reject, planned change: the plan scope adds the column, and NEWS states it.
- prior-review #1: three added lines in `R/multiverse.R` and `R/tidy.R` exceed 80 columns — reject, style: a formatter catches this, and both files were not air-clean on main (checked in T1), so the M010 lesson says not to run air on them.
- prior-review #2: `summary(direction = "yx")` has no `expect_invisible()` — fix now: add it, fixed e16657d.

After the fixes: with the R changes of e16657d stashed, the new tests fail (4 failures, 4 errors, and the tibble and min warnings). With them, `devtools::test()` gives 1984 pass, 0 fail, `devtools::document()` gives no diff, and the spell check is clean.

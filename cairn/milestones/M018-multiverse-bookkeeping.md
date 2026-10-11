# M018: Multiverse bookkeeping: repeated cells, surrogate counts, circular skips

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Resolves:** —
- **Surface tier:** user-facing — changes the grid, print, glance, and robustness counts of `synchrony_multiverse()` and the result of `autotune_wcc()`
- **Branch/PR:** m018-multiverse-bookkeeping

## Goal

Make `synchrony_multiverse()` and `autotune_wcc()` report each specification once, with the surrogate count and parameters it used, and skip circular cells on short series instead of stopping.

## Scope

**In:**
- Drop cells that repeat another cell in samples, keep the first, and announce the count. The lag cap at `floor(window / 2)` and the rounding of seconds to samples make such repeats.
- Add a per-cell `n_surrogates` grid column. `print()` says "up to" when cells differ. `glance()` gains `n_surrogates_min`.
- For wcc and wdtw on a series too short for the largest lag of the grid, skip the circular cells that cannot shift and compute the rest.
- `autotune_wcc()` records the `statistic`, `surrogate_method`, and `increment_pct` of the winning cell.
- Store `window_size`, `window_increment`, and `n_windows` as integer grid columns. `build_surface_grid()` returns an integer `n_r`.
- Roxygen and NEWS.

**Out:**
- The minimum circular shift for wgranger cells, `round(n / 4)`, stops the call or leaves few shifts for most series lengths. It goes to the new `[high]` candidate row "Granger circular surrogates in the multiverse".
- Other multiverse candidate rows (pseudo-dyad cells, IAAFT `max_iter`) stay rows.
- The unclassed "Limited unique shifts" warning of `generate_surrogate_circular()` stays as on main, also when the fallback bound of AC5 triggers it.
- `select_specification()` lines cells up across dyads by row index and does not check that the grids match. A check goes to the new `[low]` candidate row "select_specification() grid check".

## Acceptance criteria

- [ ] AC1: `synchrony_multiverse()` keeps one cell, the first in grid order, for each distinct combination of `window_size`, `lag_max`, and `window_increment` (in samples), `surrogate_method`, and `statistic`, with NA equal to NA. When it drops cells, it gives one cli message of class `bsync_duplicate_cells` that states the number dropped. `settings$n_cells`, the robustness summary, and the `n_cells_gated` of `autotune_wcc()` count only kept cells. `autotune_wcc()` gives the message once per call, not once per dyad. Tests in `test-multiverse.R` and `test-autotune.R` cover these cases: a wcc lag-cap repeat (80 Hz, `window_sec = c(1, 2)`, `lag_sec = c(0.5, 1)`, 3 cells stay, the message states 1), a window rounding repeat (`window_sec = c(1, 1.004)` at 80 Hz), an `increment_pct` rounding repeat, a wgranger repeat, a segment lag-cap repeat, a grid with no repeats and no message, and the autotune message count.
- [ ] AC2: For the pinned grids, the kept cells equal the output of main at ab7c509 for the same call and seed in `observed`, `null_mean`, `null_sd`, `es`, `p`, and `n_windows` (tolerance 1e-10). The pinned grids are a wcc grid with a lag-cap repeat and the `"phase"`, `"circular"`, and `"segment"` methods, and a wgranger grid with the `"phase"` method and no repeats.
- [ ] AC3: The grid has an integer `n_surrogates` column. It holds the number of surrogate series (matrix columns, repeated circular shifts counted) that the test of the cell used. It is `0L` wherever `n_windows` is 0. A computed cell is one with `n_surrogates > 0`. The `window_size`, `window_increment`, `n_windows`, and `n_surrogates` columns are integer for wcc, wdtw, and wgranger grids. In one call with `n_surrogates = 50`, a test asserts that a segment cell with 4 segments holds 23 and a phase cell holds 50. Tests also assert the types of these four columns for each estimator.
- [ ] AC4: `print()` of a `bsync_multiverse` shows "Surrogates per cell" as the requested count when every computed cell used it. Otherwise it shows "up to <requested> (fewest <m>)". `glance()` keeps `n_surrogates`, the requested count, and adds `n_surrogates_min`, the fewest surrogates that a computed cell used. `n_surrogates_min` is `NA` when no cell was computed or the grid has no `n_surrogates` column. For a grid without the column, `print()` shows the requested count. Tests cover both print forms, both glance cases, and the print and glance of a result whose grid lacks the column.
- [ ] AC5: For wcc and wdtw, `synchrony_multiverse()` does not stop because of the circular shift bound. A circular cell that has room for its window (n >= `window_size` + 2 * `lag_max`) but no valid shift (n <= 4 * `lag_max`) is skipped (NA values, `n_windows` 0) with a warning of class `bsync_circular_skipped`. A circular cell too short for its window is reported only as too short, as on main. The other circular cells use the largest `lag_max` among the circular cells with a valid shift as the shift bound. If n is above 4 times the largest `lag_max` of the grid, circular surrogates are made as on main (the circular pins of AC2). `autotune_wcc()` gives one `bsync_circular_skipped` warning per call, and the warning names the dyads. Tests cover three cases: n = 20 with `window_size` 10 and `lag_max` 5 (the cell is skipped and the call returns), a grid where a circular cell with a smaller lag still computes, and the autotune warning.
- [ ] AC6: The `autotune_wcc()` result has `statistic`, `surrogate_method`, and `increment_pct` fields equal to the winning grid row, with `increment_pct` NA for a segment cell. `print()` shows them. A test over a grid with two statistics and two surrogate methods asserts that the fields equal the `best_row` of `select_specification()`. A test over a segment-only grid asserts `increment_pct` is NA.
- [ ] AC7: The roxygen of `synchrony_multiverse()`, `glance.bsync_multiverse()`, and `autotune_wcc()` states the repeated-cell rule, the `n_surrogates` column, `n_surrogates_min`, the circular skip, and the new autotune fields. NEWS.md has entries for the repeated-cell drop and its effect on cell counts and `n_cells_gated`, the `n_surrogates` column, the "up to" print, `n_surrogates_min`, the circular skip, the autotune fields, and the integer column types. `devtools::test()` passes, and `devtools::check()` gives 0 errors and 0 warnings.

## Coverage

- AC1 → T3
- AC2 → T2, T3
- AC3 → T1, T4
- AC4 → T4
- AC5 → T5
- AC6 → T6
- AC7 → T3, T4, T5, T6, T7

## Tasks

- [x] T1: Integer counts. `build_surface_grid()` returns `n_r` as integer (`R/surface.R:130`, `:155`). `synchrony_multiverse()` keeps `w_samp_all` and `inc_samp_all` integer (`R/multiverse.R:329-346`, and the "Doubles" comment at `:444`). Add type tests for each estimator, and invert the `"double"` assertion at `tests/testthat/test-multiverse-segment.R:284-285`.
- [x] T2: On main at ab7c509, run the AC2 grids with a fixed seed. Pin their values in a new regression test, with a comment that names the commit and the call. The test fails until T3 drops the repeated cell. At review, re-derive the pins from ab7c509 in a separate worktree.
- [x] T3: Drop repeated cells in `synchrony_multiverse()` after the seconds-to-samples step (`R/multiverse.R:328-349`). Subset every per-cell vector together, before the segment skip warnings and surrogate generation. Use message class `bsync_duplicate_cells`. In `autotune_wcc()`, muffle the message for each dyad and give it once (`R/autotune.R:248-271`). Update the "Grid construction" roxygen.
- [x] T4: Fill the `n_surrogates` column from `ncol()` of the surrogate matrix of each cell. Update `print.bsync_multiverse()` (`R/multiverse.R:792-798`) and `glance.bsync_multiverse()` (`R/tidy.R:157-174`), with the fallback for a grid without the column. Update the roxygen `@return`.
- [x] T5: Circular fallback for wcc and wdtw (`R/multiverse.R:361-366`, `:395-398`). If n <= 4 * `max(l_samp_all)`, mark circular cells that fit their window and have n <= 4 * their lag as skipped with `bsync_circular_skipped`. Generate with the largest remaining circular lag, or make no matrix if none remain. Aggregate the warning in `autotune_wcc()` as for segment skips. Update the roxygen.
- [ ] T6: Add the `autotune_wcc()` result fields and print lines (`R/autotune.R:308-343`) and the `@return` items.
- [ ] T7: Write NEWS.md entries. Check vignettes that print multiverse grids or counts (`vignettes/choosing-parameters.Rmd`). Review each changed vdiffr snapshot and log its cause. Run `devtools::document()`, `spelling::spell_check_package()`, and `devtools::check()`.

## Work log

- 2026-10-11: created by /milestone-plan. Absorbs the candidate rows "Delivered surrogate count per multiverse cell", "Duplicate multiverse cells from the lag cap", "`autotune_wcc()` result does not record the winning cell's `surrogate_method`", "Circular surrogates at the minimum length", and "Window counts are doubles".
- 2026-10-11: question set: what to plan — the multiverse bookkeeping bundle (recommended option). How to show the per-cell surrogate count — a per-cell grid column (recommended option).
- 2026-10-11: plan gate chose dropping repeated cells over flagging them in a column because the user chose dropping at the 2026-10-07 question set. Falsified by a user need to see every requested combination in the grid.
- 2026-10-11: plan gate chose a per-cell `n_surrogates` column over print-only "up to" text because the column gives the exact count in `tidy()`. Falsified by evidence that the column breaks code that binds grids across versions.
- 2026-10-11: plan chose to skip only the circular cells that cannot shift, not every circular cell, because the rest still meet the 2 * lag_max shift rule of the generator. Falsified by a null size simulation where the kept circular cells reject above nominal.
- 2026-10-11: plan chose a new `glance()` column `n_surrogates_min` over a new meaning for `n_surrogates`, so code that reads the requested count keeps working. Falsified by a user report that the requested count misleads.
- 2026-10-11: plan widened the "Window counts are doubles" row to `window_size` and `window_increment`, which are also sample counts stored as doubles (observed in a sim_dyad run on main, 2026-10-11). `lag_max` is already integer.
- 2026-10-11: plan moved the wgranger circular shift bound to a new `[high]` candidate row. On main, `round(n / 4)` stopped the call at n = 200, 203, and 204 and left few shifts at 201 and 202 (window 2 s at 10 Hz, 19 surrogates). A new bound changes Granger circular results, so it needs its own design choice.
- 2026-10-11: inbox sweep found no open issues or outside PRs on jmgirard/bsync.
- 2026-10-11: criteria audit (full mode, fresh Opus reader) returned 10 findings, all fixed. AC5 now skips only cells that fit their window but have no valid shift, and says "does not stop because of the circular shift bound". AC2 is narrowed to the pinned grids (adds segment, wgranger uses phase), with the pin re-derivation moved to T2. AC3 defines a computed cell and sets 0L wherever `n_windows` is 0. AC1 adds `n_cells_gated`, NA keys, and window, wgranger, and segment repeat cases. AC4 adds the print fallback. AC6 adds a segment-only test. AC7 lists the NEWS entries. T1 names the test to invert. The "Limited unique shifts" warning and the `select_specification()` grid check went to Out.
- 2026-10-11: T1 done. `n_r`, `window_size`, and `window_increment` are integer. Five type-only assertions in `test-multiverse-segment.R` moved from double to integer. Suite: only the T2 pin test fails.
- 2026-10-11: T2 done. `test-multiverse-pins.R` pins 9 wcc cells (phase, circular, segment) and 2 wgranger cells from the unchanged code at ab7c509 (seeds 20261011, 20261012). The wcc pin fails until T3, because the grid still has 12 rows.
- 2026-10-11: T3 done. `inform_duplicate_cells()` gives the `bsync_duplicate_cells` message with an `n_dropped` field, and `autotune_wcc()` reads the field and gives the message once after the dyad loop. The M015 test "repeated user inputs keep their segment rows" (`window_sec = c(16, 16)`) now expects one row per specification, a planned change. Suite: 2621 pass, 0 fail, and the pins pass. No vdiffr snapshot changed.
- 2026-10-11: T4 done. The `n_surrogates` column sits after `n_windows`. Two internal helpers, `multiverse_surrogates_min()` and `multiverse_surrogate_text()`, serve `print()` and `glance()`, and both return the requested count or NA for a grid without the column. Two glance column-name tests in `test-multiverse-direction.R` gained `n_surrogates_min`. `document()` also wrote the T3 "Grid construction" text to `man/synchrony_multiverse.Rd`. Suite: 2638 pass, 0 fail.
- 2026-10-11: T5 done. `circ_skip` is computed only when "circular" is in the grid and n <= 4 * the largest lag, so other calls take the old path. With no circular cell left, the circular matrix is not made. The too-short alert count leaves out `circ_skip` cells. Tests: the n = 20 case, a mixed grid with one no-shift cell and one too-short cell (one warning), a long-series control with no warning, and the autotune warning naming dyads 1, 2, and 3. Suite: 2656 pass, 0 fail.

## Decisions

## Review

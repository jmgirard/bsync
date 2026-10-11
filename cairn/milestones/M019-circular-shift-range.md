# M019: Circular surrogates in the multiverse: a shift range per cell

- **Status:** in-progress
- **Priority:** high
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Resolves:** —
- **Surface tier:** user-facing — changes the results and warnings of the exported `synchrony_multiverse()` and `autotune_wcc()`
- **Branch/PR:** m019-circular-shift-range

## Goal

Give each circular cell of `synchrony_multiverse()` its own range of valid circular shifts, with distinct shifts only, and skip a cell that has too few shifts to reach p <= .05.

## Scope

**In:** In `synchrony_multiverse()`, a wgranger circular cell shifts `y` by at least one window (D-005). A wcc or wdtw circular cell shifts `y` by at least `2 * lag_max` for its own lag, with one surrogate matrix per distinct lag. Each cell uses `min(n_surrogates, k)` distinct shifts, where k is its number of valid shifts. A cell that fits its window and has k < 19 is skipped with a warning. This rule replaces the M018 rule (no valid shift at `4 * lag_max`). The docs, messages, and NEWS of `synchrony_multiverse()` and `autotune_wcc()` change to match. CLAUDE.md invariants 2, 6, and 7 hold: the null statistic does not change, nothing reseeds, and the result stores no shifts.

**Out:** `generate_surrogate_circular()` keeps its output and its sampling with replacement (no change planned). The wphase segment and IAAFT rows stay candidate rows. Pseudo-dyad cells in the multiverse stay a candidate row.

## Acceptance criteria

- [ ] AC1: For `estimator = "wgranger"`, `synchrony_multiverse()` draws the surrogates of each circular cell from circular shifts s with `window_size <= s <= n - window_size`. Here n is `length(y)`, and `window_size` is the window of the cell in samples. A test captures the surrogate matrix that each computed Granger circular cell receives and recovers the shift of each column by matching it to the rotations of `y`. The test uses `window_sec = c(2, 5)` at `sample_rate = 10` on series of length 200, 201, 202, 203, and 204. It asserts that range for each computed cell, and that none of the five calls stops with an error.
- [ ] AC2: For `estimator = "wcc"` and `"wdtw"`, `synchrony_multiverse()` draws the surrogates of each circular cell from circular shifts s with `2 * lag_max <= s <= n - 2 * lag_max`. Here `lag_max` is the lag of that cell in samples, not the largest lag of the grid. The test series has length 401, with `sample_rate = 1`, `window_sec = c(10, 200)`, `lag_sec = c(5, 20, 100)`, and `n_surrogates = 50`. For wcc and for wdtw, the test asserts these results. The cell with window 10 and lag 5 uses 50 distinct shifts in `10:391`, and at least one of them is outside `200:201`. The cell with window 200 and lag 20 uses 50 distinct shifts in `40:361`. The cell with window 200 and lag 100 has k = 2 and is skipped with a `bsync_circular_skipped` warning.
- [ ] AC3: A computed circular cell, for any estimator, uses `min(n_surrogates, k)` surrogates. Here k is the number of valid shifts of the cell under AC1 or AC2. No two surrogates of one cell use the same shift. The `n_surrogates` column of the grid holds that count. For wcc, wdtw, and wgranger, a test asserts the count and the distinct shifts of a cell with 19 <= k < `n_surrogates`.
- [ ] AC4: A circular cell that fits its window and has k < 19 is skipped, because no add-one p-value of its surrogates can reach p <= .05. A cell fits when n >= `window_size + 2 * lag_max`, or n >= `window_size` for wgranger. The skip also applies when the call asks for fewer than 19 surrogates. It gives a `bsync_circular_skipped` warning that names the cell and states k. The grid row of the cell has `NA` results and `n_surrogates` 0. A cell too short for its window keeps the too-short handling and gives no `bsync_circular_skipped` warning. For wcc and for wgranger, tests assert that a cell with k = 18 is skipped, that a cell with k = 19 is computed, and that a cell too short for its window gives no circular warning.
- [ ] AC5: The pinned call is `synchrony_multiverse(sim_dyad$x_A, sim_dyad$x_B, estimator = "wcc", sample_rate = 80, window_sec = c(1, 2), lag_sec = 0.5, surrogate_method = c("phase", "circular", "segment"), n_surrogates = 19)` after `set.seed(20261011)`. Its circular cells share one lag, and each has k >= 19. For every cell of that call, `observed`, `null_mean`, `null_sd`, `es`, `p`, and `n_windows` equal main at c13ffca with difference 0. `generate_surrogate_circular()` returns the same matrix and warnings as at c13ffca for the same seed in three pinned cases: `lag_max = NULL`, `lag_max = 5` on a long series, and a `lag_max` that leaves fewer valid shifts than `n_surrogates`.
- [ ] AC6: The `surrogate_method` docs of `synchrony_multiverse()` and `autotune_wcc()` state the shift ranges of AC1 and AC2, the count of AC3, and the skip of AC4, including the skip below 19 requested surrogates. The `bsync_circular_skipped` messages of both functions state the k < 19 rule. A sweep reads every hit of `grep -n -C3 -e circular -e '4 \* lag_max' -e 'n / 4'` over `R/*.R` and `vignettes/*.Rmd`. The sweep finds no doc or message text that states an old rule: one shared circular matrix, the largest lag of the grid, `n / 4` for wgranger, or a skip only at `4 * lag_max`. NEWS.md states which results change. These are circular cells of wgranger grids, of wcc and wdtw grids with more than one lag, and of cells with fewer valid shifts than `n_surrogates`. Other cells of a call with circular cells can also change, because the random draws shift. NEWS also states the new skip.

## Coverage

- AC1 → T2, T3
- AC2 → T2, T3
- AC3 → T2, T3
- AC4 → T2, T3, T4
- AC5 → T1, T4
- AC6 → T5

## Tasks

- [x] T1: In `R/surrogate_generation.R`, add an internal helper that builds circular surrogates from a minimum shift and a number of draws. Make `generate_surrogate_circular()` call it with `2 * lag_max` and keep its sampling with replacement and its warning. First, derive the three AC5 generator pins from c13ffca with `git archive` into the scratchpad (LESSONS, M018). Then add the pin tests.
- [ ] T2: In `R/multiverse.R` (about lines 423-485), give each circular cell its own shift range. Build one matrix per distinct lag (wcc, wdtw) or per distinct window (wgranger), each with `min(n_surrogates, k)` distinct shifts. Keep these draws in the circular slot of the generation order, so a single-lag call draws its random numbers as main does (AC5). Skip cells that fit and have k < 19, and state k in the warning. Remove the M018 shared-lag logic. Set `n_surrogates` per cell.
- [ ] T3: Add the tests of AC1 to AC4. Mock `.mv_wcc_cell`, `.mv_wdtw_cell`, and `.mv_granger_cell` to capture `y_surr`, and recover each shift by matching the columns to rotations of `y`.
- [ ] T4: Update the tests that relied on circular cells with few shifts. These include `test-multiverse.R` near lines 255, 605, 623, 791, 901, 926, and 953, and `test-autotune.R` near line 473. That autotune test needs a series where the lag-2 cell has k >= 19, so it still selects lag 2. Add the AC5 multiverse pin. Re-derive every pin in `test-multiverse-pins.R` from a call with circular cells, segment pins included, and record the reason in the file header. Make sure that the vdiffr snapshots do not change (their fixtures use phase).
- [ ] T5: Update the roxygen and cli messages in `R/multiverse.R` and `R/autotune.R` (about lines 85-95 and 305-320), and add the NEWS entry. Run the AC6 sweep. Run `devtools::document()`, `spelling::spell_check_package()`, `air format` on air-clean files, and `devtools::check()`.

## Work log

- 2026-10-11: created by /milestone-plan. Absorbs the candidate rows "Granger circular surrogates in the multiverse" and "Shared circular matrix with few shifts". Inbox sweep: 0 open issues, 0 open PRs.
- 2026-10-11: question set: scope is circular surrogates in the multiverse. Granger rule is at least one window. For wcc and wdtw, one matrix per lag. Cells with fewer than 19 shifts are skipped with a warning. The user approved these results changes at the gate.
- 2026-10-11: plan gate chose a minimum Granger shift of one window over 2 * ar_order, because a shift of 2 or 3 samples keeps much of the coupling in smooth signals. Falsified by a power simulation where the window rule loses power against the ar_order rule at equal size.
- 2026-10-11: plan gate chose one circular matrix per lag over one shared matrix for the largest lag, because the shared matrix leaves small-lag cells with few shifts. Falsified by evidence that cells of one grid must share their surrogate draws to be compared.
- 2026-10-11: plan gate chose to skip cells with k < 19 over computing them with a warning, because it matches the segment rule (D-002, D-004). Falsified by users who need such cells counted in the robustness summaries.
- 2026-10-11: criteria audit (full mode, fresh Opus reader) returned 9 findings, all fixed. The AC2 test now separates new code from main. Tests read shifts by capturing `y_surr`. AC1 probes two windows. D-005 records the Granger rule. AC5 is narrowed to the pinned call and three generator cases. NEWS covers changed non-circular cells. The AC6 sweep reads hits with context and covers autotune messages. AC4 keeps too-short cells apart and states the skip below 19 requested surrogates. T4 names the autotune test outcome. AC3 probes wdtw.
- 2026-10-11: lessons applied: M018 `git archive` pin derivation (T1, T4), M017 snapshots (T4), M010 air on clean files only (T5), M011 spelling (T5).
- 2026-10-11: implement started on branch m019-circular-shift-range.
- 2026-10-11: T1 done. Added internal `circular_shift_matrix()` and `draw_circular_shifts()` in `R/surrogate_generation.R`. `generate_surrogate_circular()` builds its matrix with the first. Three pins from c13ffca (seeds 101-103) pass in `test-surrogate.R`. `draw_circular_shifts()` uses `x[sample.int()]`, which is what `sample()` does for length > 1, so it draws as main does when k >= n_surrogates.

## Decisions

## Review

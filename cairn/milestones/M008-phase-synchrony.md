<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. The one size check that can fail is
     cairn_validate's <150 over the plan-owned body. -->
# M008: Phase synchrony estimator (wphase)

- **Status:** in-progress   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate -->
- **Driving RR:** —   <!-- owner: plan · create/amend-via-gate -->
- **Principles touched:** —   <!-- owner: plan · create/amend-via-gate; bsync's numbered invariants live in CLAUDE.md (1,2,4,6,7,8 bind here), no DESIGN IP/GP numbering yet -->
- **Branch/PR:** m008-phase-synchrony   <!-- owner: implement (branch) / review (PR URL) · create -->

## Goal
<!-- owner: plan · create; a wrong goal returns to plan, never edited in place -->

Add a windowed phase-synchrony estimator — `wphase()`, Hilbert instantaneous
phase via `gsignal::hilbert`, lagged phase-locking-value surface with mean
relative phase, an Rcpp core, and `wphase_surrogate()` — conforming to the
shared `bsync_surface` contract (legacy DESIGN §15 M8; promoted candidate
row). **Tier: user-facing** (exported estimator API).

## Scope
<!-- owner: plan · create/amend-via-gate -->

**In:** `src/wphase.cpp` core + `wphase()` builder (lagged surface; single
`mean_plv` aggregate statistic); NA-abort policy; `wphase_surrogate()`;
S3 methods; optima/leadership whitelist extension to `wphase`; oracle suite
per the validation doctrine incl. a frozen MNE-Python pin; `lachaux1999`
source ingestion; docs/registration. Touches `src/`: the CLAUDE.md
definition-of-done C++ clauses apply (compileAttributes, clean RcppExports,
no build artifacts), and T4 dispatches the on-demand R-hub run (recorded in
Review evidence; not a merge gate, per CLAUDE.md Git).

**Out:** wphase workflow vignette → candidate row (this commit); selectable
`statistic` for wphase (M4 pattern) → candidate row (this commit); IAAFT /
segment-shuffling / pseudo-dyad surrogates → existing M12 candidate row;
band-pass pre-filtering helpers or a narrowband-workflow guide → the
vignette candidate row carries the need; package-wide classed-condition
convention → not started here (AC5 uses house-style message matchers).

## Acceptance criteria
<!-- owner: plan · create/amend-via-gate; review reads, never reinterprets. -->

- [ ] AC1: `wphase()` is exported and returns a `bsync_surface` subclass
      (`wphase_res`) whose `results_df` carries metric columns `plv` and
      `rel_phase`; tests assert the realized window count `n_r` equals the
      hand-computed expectation under `w_max = window_size - 1`
      (Invariant 4) and the `tau` set equals a hand-enumerated expected lag
      sequence, across cases varying `window_size` (odd and even),
      `window_increment > 1`, and `lag_increment > 1`; a test supplying
      `time` asserts window positions map to real timestamps (Invariant 8);
      a test asserts the returned object's components are exactly
      `results_df` + `settings` + the aggregate — no cached phase vectors
      (Invariant 7).
- [ ] AC2: the shipped PLV/rel_phase numbers clear the validation doctrine's
      ≥2-independent-oracle-types bar via (a) closed-form — PLV and
      relative phase hand-computed from raw samples through a hand-rolled
      DFT analytic signal on one small fixed window, arithmetic in test
      comments, matched within 1e-8 using wrapped circular differences for
      phase; (c) frozen — whole-series PLV pinned within 1e-6 against
      MNE-Python on a committed fixed narrowband pair, the committed
      generator script recording package versions and reproducing bsync's
      exact preprocessing; and AC3's simulation calibration. Depth, not
      bar-clearing: (b) live — a deliberately plain pure-R windowed
      recomputation over `gsignal::hilbert` phases agreeing within 1e-8 on
      stated `sim_dyad` slices. Analytic checks: unit sinusoids at a fixed
      offset yield PLV within 1e-6 of 1 with the offset recovered in
      `rel_phase` (wrapped) within 1e-6; mean PLV on independent white-noise
      pairs stays below a committed simulated-null quantile (generator +
      seed committed). Formula source: lachaux1999, cited by the tests.
- [ ] AC3: `wphase_surrogate()` routes the observed and every surrogate
      statistic through the same aggregate helper `wphase()` uses; tests
      assert `observed_z` equals `wphase()`'s aggregate exactly and that
      passing `y` itself as the sole surrogate reproduces it (Invariant 2).
      Calibration on circular-shift surrogates, `n_surrogates = 99`, both
      blocks seeded with a same-seed-reproducibility assertion
      (Invariant 6): over 200 seeded replicates of independent white-noise
      pairs, rejections at α = .05 number between 1 and 19 inclusive
      (± 3 MC standard errors; derivation in comments); on a
      sinusoid-plus-noise pair whose stated coupling is designed for ~.95
      expected power, the rejection rate is ≥ .80. Slow parts
      `skip_on_cran`.
- [ ] AC4: `print.wphase_res()`, `summary.wphase_res()`, `plot.wphase_res()`
      (vdiffr snapshot), and `print.wphase_surr()` exist; `tidy()`,
      `glance()`, and `as_tibble()` on a `wphase_res` return the superclass
      shapes — `tidy()` the `results_df` rows; `glance()` one row whose
      columns include `window_size`, `window_increment`, `lag_max`,
      `lag_increment`, `statistic`, and the aggregate — each covered by a
      test.
- [ ] AC5: `wphase()` aborts via `cli::cli_abort()` on NA-containing input,
      asserted by a test looping over NA in `x`, NA in `y`, leading,
      interior, and all-NA cases with message matchers; the abort message
      and roxygen point at `impute_ts_gaps()` and
      `trim_edges(pad_na = FALSE)`.
- [ ] AC6: `pick_optima()` and `leadership_asymmetry()` accept the wphase
      classes (whitelists extended), each covered by a test on a `wphase`
      surface.
- [ ] AC7: every entry the NAMESPACE diff at the review ref adds —
      `export()` and `S3method()` lines — has a `_pkgdown.yml` reference
      row, and every added `export()` function carries roxygen with
      runnable `sim_dyad` examples; NEWS.md gains a plain-terms entry (no
      milestone numbers); the doc-hygiene test passes;
      `devtools::document()` produces no diff.
- [ ] AC8: at the review ref `devtools::check()` reports 0 errors,
      0 warnings, and every NOTE it emits is enumerated with a one-line
      disposition in this file's Review evidence; `devtools::test()` fully
      green.

## Coverage
<!-- owner: plan · create/amend-via-gate -->

- AC1 → T2
- AC2 → T1, T2, T3
- AC3 → T3
- AC4 → T4
- AC5 → T2
- AC6 → T4
- AC7 → T4
- AC8 → T4

## Tasks
<!-- owner: plan (create) / implement (check-off, minor edits) -->

- [x] T1: Ingest lachaux1999 (search DOI/publisher; hard stop-and-ask if
      inaccessible): source to `cairn/references/sources/`, committed
      source note + INDEX line. Write the failing oracle tests first:
      closed-form hand-DFT window fixture, sinusoid-offset and
      simulated-null-bound checks (null fixture generator + seed under
      `data-raw/`), live pure-R agreement harness. Install `gsignal`
      locally first (Imports, missing from the local library).
- [x] T2: Implement `calc_wphase_cpp` (pure core, Invariant 1; bounds → NA;
      `w_max = window_size - 1` at the boundary) and `wphase()` builder —
      `build_surface_grid(lagged = TRUE)`, settings, `mean_plv` aggregate
      helper, NA abort — greening T1; AC1's grid/time/light-object tests;
      input-validation tests; `Rcpp::compileAttributes()`; record the
      NA-policy decision as a milestone-local Decisions entry.
- [ ] T3: MNE-Python frozen pin (committed generator + fixture with
      versions and matched preprocessing); `wphase_surrogate()` sharing the
      aggregate helper; Invariant-2 identity tests; seeded circular-shift
      calibration + power tests per AC3.
- [ ] T4: Optima/leadership whitelist extension + tests; S3 methods +
      vdiffr snapshot; tidy/glance/as_tibble tests; roxygen + examples;
      `_pkgdown.yml` rows; NEWS entry; oracle-registry pointer line in
      `cairn/DESIGN.md` Conventions ("Oracle records: provenance headers in
      test files + `data-raw/` generators"); `document()` no diff; full
      `check()` clean with NOTE dispositions; dispatch the on-demand R-hub
      workflow and record its link.

## Work log
<!-- owner: any skill · append-only; one line per entry; absolute dates. -->

- 2026-08-29: created by /milestone-plan; promoted from the phase-synchrony candidate row (legacy DESIGN §15 M8).
- 2026-08-29: T1 done: lachaux1999 retrieved (Europe PMC green OA copy), source note + INDEX line committed; oracle tests authored first and run red (9 failures, wphase undefined) before any implementation; null-bound generator committed, its q99 (0.1316320, 500 reps, seed 20260829) frozen into the test.
- 2026-08-29: T3 done: MNE-Python frozen pin green (mne 1.8.0/numpy 2.0.2, generator data-raw/wphase_mne_pin.py, plv/rel_phase match within 1e-6); wphase_surrogate() + print method added sharing wphase_aggregate; Invariant-2 identity, Invariant-6 same-seed, Type-I calibration (200 reps, rejections in [1,19]) and power (design .983 pilot, floor .80) tests all green (surrogate file 83/83 with NOT_CRAN). Power design lesson: a strictly periodic common component gives the circular-shift null zero power (PLV is offset-invariant) — the test couples via shared frequency wander, recorded in test comments.
- 2026-08-29: T4 in progress: optima/leadership whitelists extended + tests; plot.wphase_res + vdiffr snapshot; tidy/glance/as_tibble tests; pkgdown rows, NEWS entry, DESCRIPTION four-estimators, DESIGN §13 oracle-records pointer; document() run (NAMESPACE + man regenerated); full suite + check() running — results at next checkpoint.
- 2026-08-29: T2 done: calc_wphase_cpp (prefix-sum cos/sin, M2 pattern) + wphase() builder + print/summary; closed-form Dirichlet oracle, sinusoid-offset, live pure-R, null-bound, grid/time/light-object, and NA-abort tests all green; full suite 1140 passing, 0 failures; RcppExports regenerated via compileAttributes.
- 2026-08-29: plan-gate criteria audit ran in full mode (user-facing tier), fresh [O] reader, two passes — pass 1: findings on every criterion (oracle-independence gap: pure-R path shared gsignal::hilbert with shipped path; missing frozen external pin vs DESIGN §13 Layer 2; unbounded tolerances; instrument-bound clauses); pass 2 over revised wording: AC1 sample-count promise unobservable → rebound to n_r + hand-enumerated lag set; white-noise bound simulated ~35% above the analytic expectation → rebound to committed simulated-null quantile; power clause a coin-flip at .80 → design-for-.95-assert-.80; glance whitelist conflict → columns named; S3method() lines escaping the NAMESPACE export() domain → domain widened; classed-condition novelty dropped for house-style message matchers; optima-layer hard-reject and src/-touching ambiguity → gate questions, both resolved (extend in scope; Rcpp core).
- 2026-08-29: plan gate chose a lagged PLV surface over lag-free windows because the surface contract and optima/leadership layer are the package's value proposition; falsified by lagged PLV proving uninterpretable or redundant with rel_phase directionality in practice.
- 2026-08-29: plan gate chose including the frozen external pin over skipping it because DESIGN §13 names Layer 2 as the only catch for misreading the field's definition; falsified by no accessible external PLV implementation with matchable preprocessing.
- 2026-08-29: plan gate chose MNE-Python over scipy+formula for the pin because a field-standard toolbox checks the convention, not just the arithmetic; falsified by MNE's PLV pipeline proving unmatchable to bsync's preprocessing on identical input.
- 2026-08-29: plan gate chose NA-abort-with-guidance over na.rm-style tolerance because the FFT-based Hilbert transform makes per-window NA handling incoherent and Invariant 3 demands loud choices; falsified by adopter demand for in-function NA handling with a coherent basis.
- 2026-08-29: plan gate chose deferring the vignette to a candidate row over parity-now to hold the 1–3 session band; falsified by users unable to adopt wphase from roxygen examples alone.
- 2026-08-29: plan gate chose extending the optima/leadership whitelists in scope over exclude-and-document because the lagged surface was chosen for that layer; falsified by phase-optima semantics proving wrong for circular data (peak PLV lag vs circular-mean issues).
- 2026-08-29: settled autonomously: Rcpp core over pure R (house architecture, DESIGN §3; Invariants 1/4/5 need the C++ boundary); falsified by the windowed-PLV loop proving cheap enough in R that the core is dead weight. Single mean_plv statistic over selectable (M4 pattern) to bound scope; falsified by immediate need for a peak statistic; the selectable form is a candidate row.

## Decisions
<!-- owner: implement / review · append-only; milestone-local -->

- 2026-08-29 (implement): NA policy — `wphase()` aborts on NA-containing
  input instead of offering `na.rm`. Rationale: the analytic signal is
  FFT-based, so one NA corrupts every phase estimate and a per-window
  `na.rm` (wcc's pattern) has no coherent meaning; Invariant 3 prefers a
  loud stop with guidance (`impute_ts_gaps()` / `trim_edges(pad_na =
  FALSE)`) over silently changing the basis. Deliberate divergence from
  wcc's `na.rm = TRUE` resolved default, settled at the M008 plan gate.

## Review
<!-- owner: review · exclusive -->

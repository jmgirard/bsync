# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Pre-migration history: see `cairn/legacy/` and git log. New IDs continue from legacy M7 (next: M016)._
_Last hygiene check: 2026-10-07 (M014 merged and archived, M011 row pruned, 3 lessons added, 3 candidate rows filed in review, validate green)_

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M015 | Segment-shuffling surrogate generator | blocked | M014 | normal | milestones/M015-segment-shuffle-surrogates.md |
| M014 | IAAFT surrogate generator | done | — | normal | milestones/archive/M014-iaaft-surrogates.md |
| M013 | Granger multiverse summary and plot for both directions | done | — | normal | milestones/archive/M013-granger-multiverse-direction.md |
| M012 | One significance rule (p <= .05) and surrogate print edges | done | — | normal | milestones/archive/M012-significance-rule.md |
<!-- rows grouped by status, not sorted by ID; keep only the 3 most recent
     terminal (done or dropped) rows — older ones live in milestones/archive/ + git -->

## Candidates

- Update CLAUDE.md for IAAFT: the "Out of scope" IAAFT bullet still calls it a candidate row, and the resolved-defaults line lists only phase and circular surrogate methods — added 2026-10-06 — M014 review (blame-history #1)
- wphase workflow vignette: parity with the wcc/wdtw/wgranger workflow articles, incl. narrowband/band-pass guidance for phase methods — added 2026-08-29 — M008 Out
- Selectable wphase aggregate statistic (the M4 wcc pattern: `statistic =` argument, null matched to observed) — added 2026-08-29 — M008 Out
- mne_connectivity PLV pipeline pin: M008's frozen pin uses MNE's analytic signal + the Lachaux formula, not mne_connectivity's own PLV pipeline (different spectral estimation; unmatchable at 1e-6); a tolerance-banded comparison against the toolbox's shipped PLV would close the convention gap the pin leaves. Promote when a definitional question about the PLV convention actually arises — added 2026-08-29 — M008 review F3
- pick_optima local search with lag_increment > 1: pick_optima_cpp hard-errors unless each window carries exactly 2*lag_max+1 lags, so `search_method = "local"` fails for any estimator surface built with lag_increment > 1 (pre-existing, all estimators; surfaced by M008's whitelist extension). Promote when a user hits it or the next optima-touching milestone — added 2026-08-29 — M008 review F13
- Enable roxygen markdown package-wide: DESCRIPTION sets no roxygen markdown option, so `[fn()]` links, `**bold**`, and backticks in existing roxygen render as raw text in `man/*.Rd`. M009's two new functions use a per-block `@md` instead. Turning it on regenerates every Rd file, so check the output. M011 added `@md` to the four `*_surrogate()` wrappers. M013 added `@md` to the four `bsync_multiverse` methods, but the `synchrony_multiverse()` block still renders its new `$robustness_yx` text raw — added 2026-10-05 — M009 T6, M013 review (diff-bug #6)
- Apply `air format` repo-wide in one formatting-only commit: `air format R tests` restyles about 26 files the change did not touch (observed 2026-10-05), so per-milestone formatting adds unrelated churn — added 2026-10-05 — M011 T4
- Wavelet coherence estimator: cross-wavelet / wavelet coherence for nonstationary time–frequency lead–lag across scales — added 2026-08-29 — legacy DESIGN §15 M9
- CRQA / MEA conventions: cross-recurrence quantification analysis and MEA-style windowed cross-correlation conveniences — added 2026-08-29 — legacy DESIGN §15 M10
- Group-level / multivariate workflow: tidy multi-dyad pipeline + aggregation / mixed-model summaries; `mvSUSY` is the >2-series reference — added 2026-08-29 — legacy DESIGN §15 M11
- Pseudo-dyad surrogates in `synchrony_multiverse()` / `autotune_wcc()`: a `surrogate_method` that draws partners from other dyads; needs the whole dyad list, which the single-dyad multiverse does not see — added 2026-10-05 — M009 Out
- CRAN submission of v0.1.0: package is submission-ready (`cran-comments.md` checklist, human-gated); release timing is user-declared, never agent-proposed — added 2026-08-29 — legacy CLAUDE.md Current focus
- Unified `bsync_ts` preprocessing object: specify before building — added 2026-08-29 — legacy DESIGN §15 unscheduled
- Expanded educational vignettes: choosing a method, interpreting a surface, reporting synchrony — added 2026-08-29 — legacy DESIGN §15 unscheduled
- [low] IAAFT on long series: at 10000 samples most surrogates reach `max_iter = 1000` before the fixed point, and `synchrony_multiverse()` / `autotune_wcc()` do not expose `max_iter`. Options: pass `max_iter` through, or stop at a spectral-accuracy target (the paper's own stop). Promote when users hit routine non-convergence warnings or slow IAAFT runs on long recordings — added 2026-10-06 — M014 review (diff-bug #2)
- [low] IAAFT rank-step tie test: the tie rule of the rank step has no test, and a test against the explicit-DFT reference is not well posed on tied input (FFT and DFT round near-ties differently). Test it through an extracted rank helper. Promote with the next change to the rank step — added 2026-10-06 — M014 review (diff-bug #8)

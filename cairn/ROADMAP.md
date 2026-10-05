# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Pre-migration history: see `cairn/legacy/` and git log. New IDs continue from legacy M7 (next: M010)._
_Last hygiene check: 2026-10-05 (no changes since 2026-08-29, validate green, no open issues, PRs, or outside merges, ROADMAP 3.3 KB, LESSONS 1.1 KB)_

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M009 | Pseudo-dyad surrogate generators | review | — | normal | milestones/M009-pseudo-dyad-generators.md |
| M008 | Phase synchrony estimator (wphase) | done | — | normal | milestones/archive/M008-phase-synchrony.md |
<!-- rows grouped by status, not sorted by ID; keep only the 5 most recent
     terminal (done or dropped) rows — older ones live in milestones/archive/ + git -->

## Candidates

- wphase workflow vignette: parity with the wcc/wdtw/wgranger workflow articles, incl. narrowband/band-pass guidance for phase methods — added 2026-08-29 — M008 Out
- Selectable wphase aggregate statistic (the M4 wcc pattern: `statistic =` argument, null matched to observed) — added 2026-08-29 — M008 Out
- mne_connectivity PLV pipeline pin: M008's frozen pin uses MNE's analytic signal + the Lachaux formula, not mne_connectivity's own PLV pipeline (different spectral estimation; unmatchable at 1e-6); a tolerance-banded comparison against the toolbox's shipped PLV would close the convention gap the pin leaves. Promote when a definitional question about the PLV convention actually arises — added 2026-08-29 — M008 review F3
- pick_optima local search with lag_increment > 1: pick_optima_cpp hard-errors unless each window carries exactly 2*lag_max+1 lags, so `search_method = "local"` fails for any estimator surface built with lag_increment > 1 (pre-existing, all estimators; surfaced by M008's whitelist extension). Promote when a user hits it or the next optima-touching milestone — added 2026-08-29 — M008 review F13
- Enable roxygen markdown package-wide: DESCRIPTION sets no roxygen markdown option, so `[fn()]` links, `**bold**`, and backticks in existing roxygen render as raw text in `man/*.Rd`. M009's two new functions use a per-block `@md` instead. Turning it on regenerates every Rd file, so check the output — added 2026-10-05 — M009 T6
- Wavelet coherence estimator: cross-wavelet / wavelet coherence for nonstationary time–frequency lead–lag across scales — added 2026-08-29 — legacy DESIGN §15 M9
- CRQA / MEA conventions: cross-recurrence quantification analysis and MEA-style windowed cross-correlation conveniences — added 2026-08-29 — legacy DESIGN §15 M10
- Group-level / multivariate workflow: tidy multi-dyad pipeline + aggregation / mixed-model summaries; `mvSUSY` is the >2-series reference — added 2026-08-29 — legacy DESIGN §15 M11
- Expanded surrogate generators: IAAFT and segment-shuffling (within-series) nulls (see `cairn/DESIGN.md` §6); the pseudo-dyad (between-dyad) generator moved to M009 — added 2026-08-29 — legacy DESIGN §15 M12 (narrowed 2026-10-05, M009 plan)
- Pseudo-dyad surrogates in `synchrony_multiverse()` / `autotune_wcc()`: a `surrogate_method` that draws partners from other dyads; needs the whole dyad list, which the single-dyad multiverse does not see — added 2026-10-05 — M009 Out
- CRAN submission of v0.1.0: package is submission-ready (`cran-comments.md` checklist, human-gated); release timing is user-declared, never agent-proposed — added 2026-08-29 — legacy CLAUDE.md Current focus
- Unified `bsync_ts` preprocessing object: specify before building — added 2026-08-29 — legacy DESIGN §15 unscheduled
- Expanded educational vignettes: choosing a method, interpreting a surface, reporting synchrony — added 2026-08-29 — legacy DESIGN §15 unscheduled

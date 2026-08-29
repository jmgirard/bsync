# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Pre-migration history: see `cairn/legacy/` and git log. New IDs continue from legacy M7 (next: M008)._
_Last hygiene check: 2026-08-29 (cairn-init migration; no milestones yet)_

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M008 | Phase synchrony estimator (wphase) | in-progress | — | normal | milestones/M008-phase-synchrony.md |
<!-- rows grouped by status, not sorted by ID; keep only the 5 most recent
     terminal (done or dropped) rows — older ones live in milestones/archive/ + git -->

## Candidates

- Phase synchrony estimator (Hilbert): analytic-signal instantaneous phase, windowed phase-locking value / phase synchrony, relative-phase output, on the shared `bsync_surface` framework — added 2026-08-29 — legacy DESIGN §15 M8. Promoted 2026-08-29 as M008 (row stays until completion)
- wphase workflow vignette: parity with the wcc/wdtw/wgranger workflow articles, incl. narrowband/band-pass guidance for phase methods — added 2026-08-29 — M008 Out
- Selectable wphase aggregate statistic (the M4 wcc pattern: `statistic =` argument, null matched to observed) — added 2026-08-29 — M008 Out
- Wavelet coherence estimator: cross-wavelet / wavelet coherence for nonstationary time–frequency lead–lag across scales — added 2026-08-29 — legacy DESIGN §15 M9
- CRQA / MEA conventions: cross-recurrence quantification analysis and MEA-style windowed cross-correlation conveniences — added 2026-08-29 — legacy DESIGN §15 M10
- Group-level / multivariate workflow: tidy multi-dyad pipeline + aggregation / mixed-model summaries; `mvSUSY` is the >2-series reference — added 2026-08-29 — legacy DESIGN §15 M11
- Expanded surrogate generators: IAAFT and segment-shuffling (within-series) plus the pseudo-dyad (between-dyad) null; pseudo-dyad depends on the multi-dyad workflow row and tests a different null (see `cairn/DESIGN.md` §2, §6) — added 2026-08-29 — legacy DESIGN §15 M12
- CRAN submission of v0.1.0: package is submission-ready (`cran-comments.md` checklist, human-gated); release timing is user-declared, never agent-proposed — added 2026-08-29 — legacy CLAUDE.md Current focus
- Unified `bsync_ts` preprocessing object: specify before building — added 2026-08-29 — legacy DESIGN §15 unscheduled
- Expanded educational vignettes: choosing a method, interpreting a surface, reporting synchrony — added 2026-08-29 — legacy DESIGN §15 unscheduled

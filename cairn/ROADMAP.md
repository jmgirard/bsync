# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Pre-migration history: see `cairn/legacy/` and git log. New IDs continue from legacy M7 (next: M012)._
_Last hygiene check: 2026-10-05 (M010 merged and archived, 1 lesson added, no new candidate rows, validate green)_

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M011 | Add-one surrogate p-values | review | — | normal | milestones/M011-add-one-p-value.md |
| M010 | Strict dyad input and named pseudo-dyad sources | done | — | normal | milestones/archive/M010-strict-dyad-input.md |
| M009 | Pseudo-dyad surrogate generators | done | — | normal | milestones/archive/M009-pseudo-dyad-generators.md |
| M008 | Phase synchrony estimator (wphase) | done | — | normal | milestones/archive/M008-phase-synchrony.md |
<!-- rows grouped by status, not sorted by ID; keep only the 5 most recent
     terminal (done or dropped) rows — older ones live in milestones/archive/ + git -->

## Candidates

- wphase workflow vignette: parity with the wcc/wdtw/wgranger workflow articles, incl. narrowband/band-pass guidance for phase methods — added 2026-08-29 — M008 Out
- Selectable wphase aggregate statistic (the M4 wcc pattern: `statistic =` argument, null matched to observed) — added 2026-08-29 — M008 Out
- mne_connectivity PLV pipeline pin: M008's frozen pin uses MNE's analytic signal + the Lachaux formula, not mne_connectivity's own PLV pipeline (different spectral estimation; unmatchable at 1e-6); a tolerance-banded comparison against the toolbox's shipped PLV would close the convention gap the pin leaves. Promote when a definitional question about the PLV convention actually arises — added 2026-08-29 — M008 review F3
- pick_optima local search with lag_increment > 1: pick_optima_cpp hard-errors unless each window carries exactly 2*lag_max+1 lags, so `search_method = "local"` fails for any estimator surface built with lag_increment > 1 (pre-existing, all estimators; surfaced by M008's whitelist extension). Promote when a user hits it or the next optima-touching milestone — added 2026-08-29 — M008 review F13
- Enable roxygen markdown package-wide: DESCRIPTION sets no roxygen markdown option, so `[fn()]` links, `**bold**`, and backticks in existing roxygen render as raw text in `man/*.Rd`. M009's two new functions use a per-block `@md` instead. Turning it on regenerates every Rd file, so check the output. M011 added `@md` to the four `*_surrogate()` wrappers — added 2026-10-05 — M009 T6
- Significance boundary for add-one p-values: print methods, `synchrony_multiverse()`, and `autotune_wcc()` call p significant at `< .05`. With (b + 1) / (n + 1), p = .05 exactly when n + 1 is a multiple of 20, for example b = 4 of 99, which b / n called significant. Decide `<` or `<=` in one place. The wphase calibration test counts `<= .05`. Promote when a user reports a boundary result or the threshold becomes an argument — added 2026-10-05 — M011 review (blame-history #1, diff-bug #12)
- Surrogate print methods at the edges: `if (x$p_value < 0.05)` crashes on an NA p-value (already on main before M011), and the four `*_surrogate()` wrappers and their print methods give no note when `n_surrogates < 20` makes p < .05 unreachable. The multiverse warns, but the wrappers do not — added 2026-10-05 — M011 review (diff-bug #2, blame-history #4)
- Apply `air format` repo-wide in one formatting-only commit: `air format R tests` restyles about 26 files the change did not touch (observed 2026-10-05), so per-milestone formatting adds unrelated churn — added 2026-10-05 — M011 T4
- Wavelet coherence estimator: cross-wavelet / wavelet coherence for nonstationary time–frequency lead–lag across scales — added 2026-08-29 — legacy DESIGN §15 M9
- CRQA / MEA conventions: cross-recurrence quantification analysis and MEA-style windowed cross-correlation conveniences — added 2026-08-29 — legacy DESIGN §15 M10
- Group-level / multivariate workflow: tidy multi-dyad pipeline + aggregation / mixed-model summaries; `mvSUSY` is the >2-series reference — added 2026-08-29 — legacy DESIGN §15 M11
- Expanded surrogate generators: IAAFT and segment-shuffling (within-series) nulls (see `cairn/DESIGN.md` §6); the pseudo-dyad (between-dyad) generator moved to M009 — added 2026-08-29 — legacy DESIGN §15 M12 (narrowed 2026-10-05, M009 plan)
- Pseudo-dyad surrogates in `synchrony_multiverse()` / `autotune_wcc()`: a `surrogate_method` that draws partners from other dyads; needs the whole dyad list, which the single-dyad multiverse does not see — added 2026-10-05 — M009 Out
- CRAN submission of v0.1.0: package is submission-ready (`cran-comments.md` checklist, human-gated); release timing is user-declared, never agent-proposed — added 2026-08-29 — legacy CLAUDE.md Current focus
- Unified `bsync_ts` preprocessing object: specify before building — added 2026-08-29 — legacy DESIGN §15 unscheduled
- Expanded educational vignettes: choosing a method, interpreting a surface, reporting synchrony — added 2026-08-29 — legacy DESIGN §15 unscheduled

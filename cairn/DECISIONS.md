# Decisions

Append-only. Never renumber; supersede with a new entry. D-entries record
choices with rationale — never deferrals ("not now" is a ROADMAP fact).

**Pre-migration decision log (pointer-only).** bsync's design decisions
through M7 live in `cairn/DESIGN.md` §14 ("Decisions resolved & remaining")
and the eight numbered Invariants in `CLAUDE.md`; the per-milestone
narrative is entombed verbatim in `cairn/legacy/MILESTONES.md`. Nothing is
re-recorded here — active work cites those anchors directly (e.g.
"DESIGN §14 #9", "Invariant 5"). New entries start at D-001.

### D-001: Surrogate p-values use the add-one form (b + 1) / (n + 1)

2026-10-05, M011 plan. Every `*_surrogate()` wrapper reports
p = (b + 1) / (n + 1), where b is the number of surrogate statistics at
least as extreme as the observed statistic and n is the number of
surrogates. This form replaces b / n. Rationale: b / n can be 0, and a test
that uses it exceeds its nominal size for most practical levels. The
add-one form is the exact Monte Carlo p-value and keeps the nominal size
(phipson2010, pp. 4, 6). The cost is that the smallest p-value is
1 / (n + 1), so a pseudo-dyad test needs at least 21 dyads to reach
p < .05. The user approved this default change at the M011 plan gate.

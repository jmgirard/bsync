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

### D-002: Surrogate significance calls use p <= .05 (narrows D-001's 21-dyad floor)

2026-10-06, M012 plan. Every significance call that bsync makes on an
add-one surrogate p-value counts p <= .05 as significant, through one
internal helper. Rationale: the size of a Monte Carlo test is
P(p <= alpha), and the add-one p-value keeps that at or below alpha
(phipson2010, pp. 4, 6). With `<`, a result whose p is exactly .05, for
example b = 4 of n = 99, is not significant, and the test is more
conservative than its nominal size. This changes the floors that D-001
states. A pseudo-dyad test with `keep_roles = TRUE` needs at least 20
dyads, not 21. The few-surrogates warning starts below 19 surrogates, not
below 20. Per-window parametric Granger p-values in `wgranger()` are not
add-one p-values and keep `< .05`. The user chose `<=` at the M012 plan
gate.

### D-003: IAAFT surrogates keep the spectrum exact by default (match = "spectrum")

2026-10-06, M014 implement stop. `generate_surrogate_iaaft()` returns the
last spectrum-adjusted series by default, so its Fourier amplitudes equal
those of `y` and its values are close. The rank-ordered series of
Schreiber & Schmitz (1996), with exact values and a close spectrum, is
`match = "values"`. Rationale: the rank step lowers the autocorrelation, so
a synchrony statistic that grows with autocorrelation reads high against
rank-ordered surrogates, and a WCC test against them rejected a true null
too often in the M014 size simulation. The spectrum form kept the nominal
size. The measured rates are in the M014 milestone record. The
`"iaaft"` method of `synchrony_multiverse()` and `autotune_wcc()` uses
this default. The user chose it at the M014 implement stop.

### D-004: Segment surrogates use uniform non-identity orders and aligned windows

2026-10-07, M015 RR01 ingest. `generate_surrogate_segment()` draws
distinct segment orders uniformly from all orders except the original, so
a surrogate can leave a segment in place. A segment test is run only with
windows where every lagged window lies inside one segment: segment size at
least window + 2 * lag_max with the increment equal to the segment size,
or segment = window = increment for `wgranger()`. The multiverse forces
this design for `"segment"` cells. Rationale: the add-one p-value is exact
only when the surrogate orders and the identity form a group and the
segments of `y` are exchangeable (phipson2010, p. 6). Orders with no
segment in place are not such a group, and windows that cross segment
boundaries keep continuity in the observed series that the surrogates
lose. Both made the test reject a true null too often. The measured rates
are in RR01 (`cairn/reviews/archive/RR01-segment-surrogate-size.md`). This
replaces the M015 plan choice of orders with no segment in place and of
segment size equal to the window size.

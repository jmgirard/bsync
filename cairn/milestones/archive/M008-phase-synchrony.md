# M008: Phase synchrony estimator (wphase)

**Status:** done (2026-08-29, PR #3 https://github.com/jmgirard/bsync/pull/3)

**Goal:** Add a windowed phase-synchrony estimator — `wphase()` (Hilbert
phase, lagged PLV surface + mean relative phase, Rcpp core) and
`wphase_surrogate()` — on the shared `bsync_surface` contract.

**Outcome:** `wphase()` + `calc_wphase_cpp` (prefix-sum cos/sin core),
`wphase_surrogate()` (shared aggregate helper, Invariant 2), NA-abort policy,
`print/summary/plot.wphase_res`, `print.wphase_surr`,
`print/summary.wphase_optima`, optima/leadership whitelists extended,
tidy/glance/as_tibble coverage. Oracles: closed-form Dirichlet fixture,
frozen MNE-analytic-signal pin (data-raw/wphase_mne_pin.py), live pure-R,
simulated-null bound (data-raw/wphase_null_bound.R), Type-I calibration +
power tests. lachaux1999 ingested (references/lachaux1999.md); DESIGN §4
four-estimator table + §13 oracle-records pointer; NEWS/README/pkgdown.

**Decisions:** NA policy — abort with guidance instead of `na.rm` (FFT-based analytic signal; Invariant 3); milestone-local, in git history.

**Review:** three-lens fan-out. Blame: none; prior-review: 1 (WORDLIST,
fixed); diff-bug [O]: core brute-force-verified to 1.2e-15, 17 findings —
14 fixed at the gate (backwards broadband caveat, saturated snapshot,
overlay mislabel, missing optima methods among them), 2 routed as candidate
rows, 2 rejected with reasons. check() 0/0/0; CI 8/8; suite 1172 passing.

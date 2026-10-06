# M011: Add-one surrogate p-values

**Status:** done (2026-10-06, PR #6 https://github.com/jmgirard/bsync/pull/6)

**Goal:** Make every surrogate test report the add-one p-value (b + 1) / (n + 1)
(phipson2010 (p. 6)), so that a reported p-value is never 0 and the test keeps
its nominal size.

**Outcome:** The p-value lines in `wcc_surrogate()`, `wdtw_surrogate()`,
`wgranger_surrogate()` (both directions), and `wphase_surrogate()` use
(b + 1) / (n + 1). The wrappers abort on a zero-column surrogate matrix.
Internal `format_p_value()` prints 4 digits in fixed notation or `< 0.0001`.
`warn_few_surrogates()` (class `bsync_few_surrogates`) warns below 20 in
`synchrony_multiverse()`, and `autotune_wcc()` warns once and muffles the
per-dyad copies. Source note `cairn/references/phipson2010.md`, DESIGN §6 and §9,
`@md` on the four wrapper docs, vignettes, README, and NEWS are updated.

**Decisions:** D-001 (add-one form). The plan-gate approval stands as the
pre-1.0 waiver of a deprecation cycle.

**Review:** No returns. Claim audit: 86 claims, 4 corrected. Three-lens
fan-out: 22 findings, with 12 fixed and 7 rejected. 3 went to candidate rows
("Significance boundary for add-one p-values", "Surrogate print methods at
the edges"). check() 0/0/0, CI 9/9.

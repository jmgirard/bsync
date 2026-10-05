# M009: Pseudo-dyad surrogate generators

**Status:** done (2026-10-05, PR #4 https://github.com/jmgirard/bsync/pull/4)

**Goal:** Add two exported generators that build pseudo-dyads (one partner's
series paired with a partner from a different dyad), so users can test
synchrony against the between-dyad pseudo-synchrony null.

**Outcome:** `generate_surrogate_pseudo()` gives a per-dyad partner matrix for
all four `*_surrogate()` wrappers. It excludes shorter partners with a warning
and crops longer ones start-aligned. `generate_pseudo_dyads()` gives the
sample-wide N(N-1) or 2N(N-1) pairings as `list(x, y)`. Both carry a
`"sources"` attribute and default to `keep_roles = TRUE`. Source note:
references/kleinbub2020.md (rMEA 1.2.2 `shuffle()`). Also a vignette section
with a worked example, NEWS, README reference, pkgdown rows, DESIGN §2/§6.

**Decisions:** none cross-cutting. The plan-gate choices (start-aligned crop,
roles kept by default) are in the milestone file's work log in git history.

**Review:** two returns (missing `@seealso` links, plan body at the 150-line
cap), both fixed. Claim audit: 66 claims, 3 corrected. Three-lens fan-out:
22 findings, 13 fixed, 3 to candidate rows (strict `.extract_xy()`, dyad
names in `sources`, add-one p-value), 4 rejected, 2 noted. check() 0/0/0,
CI 9/9 green.

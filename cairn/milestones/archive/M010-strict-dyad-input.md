# M010: Strict dyad input and named pseudo-dyad sources

**Status:** done (2026-10-05, PR #5 https://github.com/jmgirard/bsync/pull/5)

**Goal:** Make every function that reads a `dyad_list` identify each dyad's
series by exact name or by position, never by partial match or a silent
role swap, with errors that name the dyad, and carry `dyad_list` names into
the pseudo-dyad generators' output.

**Outcome:** `.extract_xy(dyad, index)` reads `x`/`y` by exact name when
each occurs once, else the first two elements, and aborts on any other x/y
naming with the dyad index. `autotune_wcc()` extracts every dyad before
sampling, keeps `dyad_list` names on `$dyad_multiverses`, and rejects a
single data frame. `.pseudo_dyad_names()` validates names (non-empty, not NA,
unique). Both generators add `name` / `x_name` / `y_name` columns to
`"sources"`, plus column names and element names when named. NEWS lists the
inputs whose results change without an error. Testthat failure artifacts are
now in `.gitignore`.

**Decisions:** none cross-cutting (plan-gate choices: work log, git).

**Review:** one return (a missing AC3 probe for unnamed lists under an
integer draw), fixed. Claim audit: 44 claims, 4 corrected. Three-lens
fan-out: 4 fixed (names lost in `autotune_wcc()` output, separator labels,
two input-shape error messages), 7 rejected, 2 noted. check() 0/0/0, CI 9/9.

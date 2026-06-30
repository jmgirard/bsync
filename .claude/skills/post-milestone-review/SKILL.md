---
name: post-milestone-review
description: Read-only conformance audit of a completed milestone against CLAUDE.md acceptance criteria/invariants and the DESIGN.md contracts. Use after finishing a milestone, before starting the next.
allowed-tools: Read, Grep, Glob, Bash
argument-hint: "[milestone-number]"
---

0. Confirm $ARGUMENTS appears in CLAUDE.md's "## Completed milestones" section (e.g., look for
   "M$ARGUMENTS (done)"). Do NOT rely on "## Current focus" — a finalized milestone should have been
   moved to "Completed milestones" by /implement-milestone's finalize step. If it is NOT there, stop
   and flag the *actual* state rather than guessing:
   - Still listed as active/next in "## Current focus" **with** implementation commits on `main` →
     it was implemented but never finalized (the finalize bookkeeping is missing). Say so and point
     the user to finalize it first — this is the /implement-milestone finalize step, **not** a
     re-implementation.
   - Absent everywhere **and** no implementation commits → not implemented yet; run
     /implement-milestone first.
   - Present but numbered differently than the work actually on `main` → the milestone number is
     wrong.

# Post-Milestone Review: Milestone $ARGUMENTS

This is a **read-only audit**. Do not edit, create, or fix any files — report only.

1. Re-read `CLAUDE.md` and `DESIGN.md` fresh; don't rely on conversation memory of what they say.
2. Run the test suite and `devtools::check()` (use `--as-cran` for any release-track milestone, e.g.
   the M3 CRAN-readiness work) and report actual results, not assumptions.
3. For each acceptance criterion for Milestone $ARGUMENTS in `CLAUDE.md`: state **Met / Partially Met
   / Not Met**, with the specific file/test as evidence.
4. For each of the eight invariants in `CLAUDE.md`: state whether the implementation upholds it,
   citing the file/location. Pay special attention to:
   - **Invariant 1** — no `*_cpp` function exposed as user API; all validation in the R wrapper.
   - **Invariant 2** — the surrogate null statistic exactly matches the observed statistic.
   - **Invariant 5** — any C++ change carries a numerical-regression test; confirm it actually
     compares against a reference and would catch a divergence.
   Flag silent violations even if no test currently catches them.
5. Compiled-code hygiene specific to this package:
   - No build artifacts (`src/*.o`, `src/*.so`, `src/*.dll`, `Rplots.pdf`, `.DS_Store`) are tracked
     (`git ls-files`).
   - `RcppExports.cpp`/`RcppExports.R` are clean against a fresh `Rcpp::compileAttributes()` (no
     uncommitted reordering churn).
   - `src/Makevars*` does not link OpenMP unless `#pragma omp` is actually used (the M2 decision).
6. List untested edge cases implied by `DESIGN.md`/`CLAUDE.md` that the current suite does not
   exercise (e.g., short-series abort path, `na.rm = FALSE`, realized window length, NA-laden
   surrogate windows, coupled-vs-independent surrogate p-value sanity, Granger degrees-of-freedom
   guard).
7. Dependency hygiene: anything newly in `Imports` that should be `Suggests` (or vice versa); any
   new dependency added without a flag; confirm heavy compute stayed in the C++ cores.
8. Roxygen check: does documentation explain *why* defaults were chosen (per CLAUDE.md's
   documentation standard), not just *what* — including the `suggest_wcc_params()` heuristic and the
   surrogate-null semantics.
9. Output a single triaged list: **Blocking** (must fix before the next milestone) / **Should-fix** /
   **Nice-to-have**. End with one line: **READY** or **NOT READY** for the next milestone, and why.

Do not rewrite, fix, or refactor anything during this review — that's a separate step.

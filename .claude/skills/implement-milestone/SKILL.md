---
name: implement-milestone
description: Implements the milestone plan already established in this conversation, following CLAUDE.md's dev workflow and definition of done. Use right after planning a milestone in this same session.
disable-model-invocation: true
argument-hint: "[milestone-number]"
---

0. Confirm $ARGUMENTS matches the milestone number recorded as active in CLAUDE.md's "## Current
   focus" section. If it doesn't match, stop and flag the mismatch rather than proceeding.

# Implement Milestone $ARGUMENTS

Implement Milestone $ARGUMENTS using the plan already established earlier in this conversation.

**If no milestone plan is visible in the current context, stop and ask for it — do not invent or
assume one.**

Follow `CLAUDE.md`'s dev workflow and definition of done throughout:
- Build in small, reviewable increments: one coherent change → test → commit → repeat. Do not produce
  one large commit at the end.
- **After any C++ (`src/`) edit:** recompile via `devtools::load_all()`; if a `[[Rcpp::export]]`
  signature changed, run `Rcpp::compileAttributes()` and commit the regenerated
  `RcppExports.cpp`/`RcppExports.R` **clean** (no stray reordering churn). Never stage or commit build
  artifacts (`*.o`, `*.so`, `*.dll`).
- **Any C++ optimization** must ship with a numerical-regression test asserting output matches the
  prior implementation within tolerance on `sim_dyad` and an NA-laden case (Invariant 5). Performance
  must not change results.
- **Efficiency milestones:** record before/after timings with the `bench/` script and state the
  measured speedup; "it's faster" is not acceptance.
- After any roxygen change, run `devtools::document()` and commit the regenerated `NAMESPACE`/docs.
- Write or update tests for new behavior; run `devtools::test()` before each commit.
- Run `devtools::check()` (`--as-cran` on release-track milestones) before considering any sub-step
  done; resolve errors/warnings, triage notes.
- Style and lint (`air format` or `styler::style_pkg()`, `lintr::lint_package()`).
- Respect the ask-first guardrails in `CLAUDE.md`: flag before adding an `Imports` dependency,
  changing a resolved default, changing the numeric output of a C++ core (beyond a deliberately
  planned semantics fix), adopting OpenMP, or touching git history/tags.
- If you hit a design ambiguity not covered by `DESIGN.md`, stop and ask rather than deciding
  unilaterally — design decisions belong in `DESIGN.md`.
- Update `NEWS.md` for user-visible changes (create it if it does not yet exist).

Do not run `/post-milestone-review` yourself at the end — that's a separate, deliberate step the user
triggers after reviewing the work.

---
name: implement-milestone
description: Implements the milestone plan already established in this conversation, following CLAUDE.md's dev workflow and definition of done. Use right after planning a milestone in this same session.
disable-model-invocation: true
argument-hint: "[milestone-number]"
---

0. Confirm $ARGUMENTS matches the milestone number recorded as active in CLAUDE.md's "## Current
   focus" section. If it doesn't match, stop and flag the mismatch rather than proceeding.
0b. Confirm you are on the milestone's feature branch (created by /plan-milestone, named
    `m$ARGUMENTS-<slug>`), not on `main`. If you are on `main`, stop and flag — milestone work must
    not be committed directly to `main`.

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

## Finalize (only when the milestone is genuinely done)

When — and only when — **every** acceptance criterion for this milestone in `CLAUDE.md` is
demonstrably met **and** `devtools::check()` (`--as-cran` on release-track milestones) is clean
(0 errors / 0 warnings, notes triaged) **and** the working tree is committed, finalize the
bookkeeping yourself (this is mechanical, not a judgment call — the human quality gate is
`/post-milestone-review`, run separately). **`MILESTONES.md` is the single source of truth for
milestone history; do not duplicate the narrative across files:**

- Add a detailed `M$ARGUMENTS (done)` entry to **`MILESTONES.md`**, placed **in numeric order**
  (immediately after `M{$ARGUMENTS − 1}`; never append out of order or skip a number), in the same
  `(done)` style as the existing entries: the met acceptance criteria, the commit range, the test
  count, and the `R CMD check` result; note any post-plan deviations.
- Add the matching **one-line index entry** to CLAUDE.md's "## Completed milestones" list (numeric
  order), and clear/replace "## Current focus" so the next milestone is at the front (mark it "next").
- Record any DESIGN.md **contract** changes in the relevant numbered section (§1–14) — **not** as a
  second copy of the milestone log (DESIGN.md §15 is only the forward roadmap + a pointer to
  `MILESTONES.md`).
- Do **not** duplicate the milestone narrative across files: detail in `MILESTONES.md`, one line in
  CLAUDE.md, user-facing notes in `NEWS.md`, design-contract deltas in DESIGN.md §1–14.
- Commit this bookkeeping **on its own** (on the feature branch), not bundled with any implementation
  commit (e.g. `Mark M$ARGUMENTS complete; log to MILESTONES.md, index in CLAUDE.md`).

If any acceptance criterion is unmet, the check is not clean, or a guardrail decision is still
pending, **stop and report** — do **not** mark the milestone complete or write the `MILESTONES.md`
entry.

## Finishing the milestone (branch → PR, not direct to `main`)

- All milestone commits live on the `m$ARGUMENTS-<slug>` feature branch. **Push the branch and open
  a PR into `main` only when the user asks** (respects CLAUDE.md's "push only when asked" cadence) —
  never auto-push.
- The PR merges into `main` only after the three GitHub Actions workflows go green: `R-CMD-check`,
  `test-coverage`, `pkgdown`. Do not merge on red CI.
- **If this milestone touched `src/` (a C++ core), recommend an on-demand R-hub run**
  (`rhub::rhub_check()`, workflow `.github/workflows/rhub.yaml`) for sanitizers/valgrind/extra
  platforms before merge — it is not part of the per-PR gate but is the highest-value check for the
  Armadillo cores. Likewise flag it before any actual CRAN submission.
- Trivial, isolated doc-typo fixes may still go directly to `main` at the user's discretion; anything
  touching `R/`, `src/`, `tests/`, `DESCRIPTION`, or vignettes goes through the PR.

Do not run `/post-milestone-review` yourself at the end — that's a separate, deliberate step the user
triggers after reviewing the work.

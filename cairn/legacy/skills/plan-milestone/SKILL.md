---
name: plan-milestone
description: Scaffolds a milestone plan against CLAUDE.md and DESIGN.md, prompting for milestone-specific details before producing the plan.
disable-model-invocation: true
argument-hint: "[milestone-number]"
---

# Plan Milestone $ARGUMENTS

1. Read CLAUDE.md's "## Current focus" section to confirm the currently/previously recorded
   milestone number. If $ARGUMENTS isn't exactly one more than the last completed milestone (or
   doesn't otherwise make sense — e.g., re-planning the same milestone), stop and flag the
   discrepancy before proceeding.

2. Ask the user for milestone-specific details: what this milestone covers, why it comes next, and
   any constraints, deviations from DESIGN.md, or open questions they already have in mind. Wait for
   their answer — do not guess or proceed on assumptions.

Once they respond:

3. Re-read CLAUDE.md and DESIGN.md fresh. For prior-milestone detail, the single log is
   `MILESTONES.md` (CLAUDE.md's "## Completed milestones" carries only a one-line index); the forward
   roadmap is DESIGN.md §15 and deferred/remaining decisions are DESIGN.md §14.
4. Confirm this milestone's scope against DESIGN.md §15 (roadmap) and §14 (decisions remaining), and
   CLAUDE.md's "out of scope for now" list; confirm it isn't already done per the `MILESTONES.md`
   log; flag if the user's details conflict with any of these.
5. Propose a concrete plan: files to create/modify (name the R wrappers, C++ cores, tests, and docs),
   the order of implementation, and testable acceptance criteria in the same style as M1's (see
   `MILESTONES.md`).
6. For any milestone touching a C++ core, the plan must include: a `Rcpp::compileAttributes()`
   regeneration step, a **numerical-regression test** vs. a reference implementation on `sim_dyad`
   (Invariant 5), and — for efficiency work — a `bench/` before/after timing step with a measured
   target.
7. Flag any design ambiguity the brief doesn't resolve — do not invent a resolution; if it's a design
   decision, it belongs in DESIGN.md, so ask.
8. Do not write any code. This is planning only; implementation happens via /implement-milestone
   afterward.

9. Once the user explicitly approves the plan (branch → PR workflow; milestone work never lands
   directly on `main`):
   a. Create the milestone's feature branch off an up-to-date `main`
      (e.g. `git switch main && git pull && git switch -c m$ARGUMENTS-<short-slug>`). All milestone
      work — the planning commit and implementation — lands on this branch; `main` is updated only
      by merging a reviewed PR once CI is green (see /implement-milestone).
   b. On the branch, update CLAUDE.md's "## Current focus" section to make Milestone $ARGUMENTS active
      (with its acceptance criteria), and commit that change **on its own** (not bundled with any
      implementation commit). Commit **locally only** — do not push the branch or open the PR until
      the user asks (respects CLAUDE.md's "push only when asked" cadence).

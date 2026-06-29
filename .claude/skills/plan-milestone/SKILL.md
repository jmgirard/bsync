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

3. Re-read CLAUDE.md and DESIGN.md fresh.
4. Confirm this milestone's scope against DESIGN.md §15 (roadmap) and §14 (decisions remaining), and
   CLAUDE.md's "out of scope for now" list; flag if the user's details conflict with either.
5. Propose a concrete plan: files to create/modify (name the R wrappers, C++ cores, tests, and docs),
   the order of implementation, and testable acceptance criteria in the same style as M1's in
   CLAUDE.md.
6. For any milestone touching a C++ core, the plan must include: a `Rcpp::compileAttributes()`
   regeneration step, a **numerical-regression test** vs. a reference implementation on `sim_dyad`
   (Invariant 5), and — for efficiency work — a `bench/` before/after timing step with a measured
   target.
7. Flag any design ambiguity the brief doesn't resolve — do not invent a resolution; if it's a design
   decision, it belongs in DESIGN.md, so ask.
8. Do not write any code. This is planning only; implementation happens via /implement-milestone
   afterward.

9. Once the user explicitly approves the plan, update CLAUDE.md's "## Current focus" section to make
   Milestone $ARGUMENTS active (with its acceptance criteria), and commit that change on its own (not
   bundled with any implementation commit).

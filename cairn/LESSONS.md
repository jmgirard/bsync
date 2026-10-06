# Lessons

Durable repo lessons — build quirks, testing tricks, gotchas worth
remembering next time — captured at milestone end and surfaced at plan time.
Not status, not decisions: a lesson is a reusable "how this repo actually
behaves" note. Cross-cutting *choices* still go to `DECISIONS.md`.

One line per lesson: `- YYYY-MM-DD (M<NNN>): <lesson>`. Two caps: 50 lines
and 20,000 bytes; over either, retire or prune before adding. Corrected in
place when proven false (never append a correction).

- 2026-08-29 (M008): PLV is invariant to constant phase offsets, so a circular-shift null has ZERO power against a strictly periodic shared component — power/coupling test signals must carry shared frequency wander, not a pure sinusoid.
- 2026-08-29 (M008): cli output (cli_h1/cli_dl/alerts) lands on the message stream under testthat — assert print methods with expect_message or snapshots, never expect_output.
- 2026-10-05 (M009): roxygen markdown is off package-wide (no DESCRIPTION option), so `[fn()]` links, `**bold**`, and backticks render raw in man/*.Rd. Add `@md` to new roxygen blocks until the package-wide candidate row lands.
- 2026-10-05 (M010): `R/autotune.R` and `tests/testthat/test-autotune.R` are not air-clean on main, so `air format` on them restyles untouched code. Run `air format` only on files that were air-clean before the edit, and stage explicit paths, not `git add -A`.
- 2026-08-29 (M008): a vdiffr fixture over a saturated surface (all values ~1) is blind to fill-mapping changes — build snapshot fixtures with real dynamic range (e.g. different-frequency signal pairs for PLV).

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
- 2026-08-29 (M008): a vdiffr fixture over a saturated surface (all values ~1) is blind to fill-mapping changes — build snapshot fixtures with real dynamic range (e.g. different-frequency signal pairs for PLV).

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
- 2026-10-06 (M011): `future` multisession workers load the installed bsync, not the `load_all()` copy, so vignette code run against dev code fails with "no package called 'bsync'". Swap in `plan(sequential)`. Results match, because surrogate draws use `future.seed = TRUE`.
- 2026-10-06 (M011): the spelling test compares against `tests/spelling.Rout.save`, so a new proper noun in roxygen, NEWS, or a vignette gives a check NOTE until it is added to `inst/WORDLIST`. Run `spelling::spell_check_package()` before `devtools::check()`.
- 2026-10-06 (M012): cli alert functions read their first argument as a glue template, so a variable message that holds `{` errors ("object 'b' not found"). Pass variable text as `cli::cli_alert_info("{msg}")`, never as the template itself.
- 2026-10-07 (M013): `vdiffr::expect_doppelganger()` skips under `testthat::test_file()` from Rscript, because `NOT_CRAN` is unset. Record or update a snapshot with `devtools::test(filter = "<file>")`, which sets `NOT_CRAN=true`.
- 2026-10-07 (M014): a surrogate generator can pass every invariant test and still give a liberal test. The rank-ordered IAAFT form kept exact values and a near spectrum, yet rejected 10.5% of null pairs at a nominal 5%. Check each new generator with a null size simulation that uses a real synchrony statistic.
- 2026-10-07 (M014): the cairn merge guard reads `gh pr merge <N>` only as the whole command. A `cd` prefix in another form, or a second command joined after it, is denied. Run the merge alone from the repo directory.
- 2026-10-07 (M014): cli `{?s}` and `{?has/have}` after a numeric vector pluralize on its value, not its length, so `Dyad{?s} {idx}` with `idx = 2` prints "Dyads". Write `{cli::qty(length(idx))}` before the pluralized word.
- 2026-10-07 (M015): a permutation surrogate set gives an exact add-one p-value only if its orders, with the original, form a group drawn uniformly. Orders that move every segment break this and rejected a true null 9-12% of the time at 5 segments. Draw from all non-identity orders, and keep every lagged window inside one segment (segment = window + 2 * lag, increment = segment).
- 2026-10-07 (M015): to plant a defect in an exported function under `devtools::load_all()`, assign the changed copy in both the namespace and `as.environment("package:bsync")`. Tests call the package-environment copy, so a namespace-only swap leaves them green. Count `r$error` as well as `r$failed` from `as.data.frame(test_file())`.
- 2026-10-07 (M015): a GitHub Actions job cancelled after hours inside `r-lib/actions/setup-r` never ran its tests. It is a runner hang, not a code failure. Re-run the job before treating the check as red.
- 2026-08-29 (M008): a vdiffr fixture over a saturated surface (all values ~1) is blind to fill-mapping changes — build snapshot fixtures with real dynamic range (e.g. different-frequency signal pairs for PLV).

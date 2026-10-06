# M012: One significance rule (p <= .05) and surrogate print edges

**Status:** done (2026-10-06, PR #7 https://github.com/jmgirard/bsync/pull/7)

**Goal:** Make every surrogate significance call in bsync come from one
internal p <= .05 rule that also handles NA p-values and unreachable thresholds.

**Outcome:** Two internal helpers set the rule: `is_significant()`
(`!is.na(p) & p <= 0.05`) and `significance_reachable()`. The four surrogate
print methods use them through `print_significance_call()` and
`print_surrogate_count_notes()`. So do `multiverse_robustness()`,
`multiverse_plot_data()`, and the `autotune_wcc()` significance rate. An NA
p-value prints "No significance call", and 18 or fewer surrogates print a note
that no result can reach p <= .05. `warn_few_surrogates()` moves from 20 to 19.
Roxygen, five `man/` pages, two vignettes, and NEWS state p <= .05, 19
surrogates, and 20 dyads. Per-window Granger F-test p-values keep `< 0.05`.

**Decisions:** D-002 (p <= .05, narrows D-001's dyad floor).

**Review:** One defect return (AC5: two wgranger NA cases lacked
`expect_invisible()`). Three-lens fan-out: 17 findings, 4 fixed in c1eceb5,
1 to a candidate row, 12 rejected. The Codecov PR comment went to a second
row ("Test the all-NA error in `plot.bsync_multiverse()`"). check() 0/0/0, CI
9/9. The last 3 work-log commits (CI-wait stop, resume 2, approval) were not
pushed and survive only in the local reflog at 2eccc2f.

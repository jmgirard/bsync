# M013: Granger multiverse summary and plot for both directions

**Status:** done (2026-10-07, PR #8 https://github.com/jmgirard/bsync/pull/8)

**Goal:** Let users read the y -> x direction of a Granger `synchrony_multiverse()` result through the print, summary, glance, and plot methods that before reported only x -> y.

**Outcome:** Granger results carry `$robustness_yx`, computed by `multiverse_robustness()` from `es_yx`/`p_yx`. The internal helper `multiverse_direction()` matches `direction = c("xy", "yx")`. It stops for `"yx"` on WCC, WDTW, or a grid without `_yx` columns. It recomputes `$robustness_yx` from the grid for older Granger objects. `print()`, `summary()`, `glance()`, and `plot()` take `direction` (last named argument in `plot()`). For Granger only, they show a "Direction" line, a `glance()` `direction` column, and a title that names the direction. The helper `multiverse_plot_title()` builds the title. `summary()` prints "ES range: none" when no ES is computable. The vignette, a header comment, and a test title no longer name a nonexistent `es_xy` column. New tests: `test-multiverse-direction.R` and the `multiverse-granger-yx` vdiffr snapshot. The "Test the all-NA error" candidate row was absorbed.

**Decisions:** none. The plan gate chose a `direction` argument over always showing both directions.

**Review:** Criteria audit at plan: 13 findings, all fixed. Claim audit: 60 claims, 4 corrected. Three-lens fan-out: 23 findings. 9 were fixed in e16657d, among them older Granger objects that lacked `$robustness_yx` and the `[Inf, -Inf]` ES range. 1 went to the candidate row "Enable roxygen markdown package-wide". 5 were rejected and 8 noted. No returns or amendments. Local check() at e6de332 gave 0/0/0. CI was 9/9 green. No lesson was retired.

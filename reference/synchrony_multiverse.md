# Synchrony Multiverse Analysis

Sweeps a seconds-specified parameter grid across analytic choices and
evaluates each specification with a matched-null surrogate test
(Invariant 2). The headline metric is \*\*effect size vs. the null\*\* –
not raw synchrony, which autocorrelation inflates.

## Usage

``` r
synchrony_multiverse(
  x,
  y,
  estimator = c("wcc", "wdtw", "wgranger"),
  sample_rate,
  window_sec,
  lag_sec = NULL,
  increment_pct = 0.1,
  statistic = "mean_abs_z",
  surrogate_method = "phase",
  n_surrogates = 100L,
  ar_order = 1L,
  scale_method = c("global", "local", "none"),
  distance_metric = c("L2", "L1"),
  na.rm = TRUE
)
```

## Arguments

- x:

  Numeric vector; the reference time series.

- y:

  Numeric vector; the query time series (same length as \`x\`).

- estimator:

  Character string; windowed estimator: \`"wcc"\` (default), \`"wdtw"\`,
  or \`"wgranger"\`.

- sample_rate:

  Single positive number; sampling rate in Hz, used to convert
  \`window_sec\` and \`lag_sec\` to samples.

- window_sec:

  Numeric vector; window size(s) in seconds to sweep.

- lag_sec:

  Numeric vector; max lag(s) in seconds. Required for \`"wcc"\` and
  \`"wdtw"\`; ignored for \`"wgranger"\`. Capped at \`window_sec / 2\`
  per cell.

- increment_pct:

  Numeric vector; window increment(s) as a fraction of \`window_size\`
  (e.g., \`0.1\` = 10% step). Default is \`0.1\`.

- statistic:

  Character vector; aggregate statistic for \`"wcc"\` only. One or both
  of \`"mean_abs_z"\` (SUSY) and \`"peak"\` (rMEA/Boker). Ignored for
  other estimators.

- surrogate_method:

  Character vector; surrogate generator(s): `"phase"` (preserves power
  spectrum), `"circular"` (preserves autocorrelation), and/or `"iaaft"`
  (preserves the power spectrum exactly and the value distribution
  closely;
  [`generate_surrogate_iaaft`](https://jmgirard.github.io/bsync/reference/generate_surrogate_iaaft.md)
  with its default `match = "spectrum"`). `"iaaft"` needs `y` without
  missing values.

- n_surrogates:

  Single positive integer; number of surrogates per cell. Default is
  \`100\`. Use \>= 1000 for reporting. A cell is significant at p \<=
  .05. Below 19 the call warns, because the smallest possible p-value, 1
  / (n_surrogates + 1), is then above .05 and no cell can be
  significant.

- ar_order:

  Single positive integer; AR order for \`"wgranger"\`. Default is
  \`1L\`.

- scale_method:

  Character string; scaling for \`"wdtw"\`. Default is \`"global"\`.

- distance_metric:

  Character string; distance metric for \`"wdtw"\`. Default is \`"L2"\`.

- na.rm:

  Logical; passed to \`"wcc"\`. Default is \`TRUE\`.

## Value

A \`bsync_multiverse\` object with:

- \`\$grid\`:

  \[tibble::tibble()\] with one row per parameter cell: specification
  columns, \`window_size\`/\`lag_max\`/\`window_increment\` (samples),
  \`n_windows\`, \`observed\`, \`null_mean\`, \`null_sd\`, \`es\`, \`p\`
  (x -\> y for Granger, plus \`observed_yx\`, \`null_mean_yx\`,
  \`null_sd_yx\`, \`es_yx\`, and \`p_yx\` for y -\> x).

- \`\$settings\`:

  Named list of call-level inputs.

- \`\$robustness\`:

  Named list: \`n_cells\` (total specifications in the grid),
  \`n_valid\` (cells with a computable ES in that direction; the rest
  were skipped as too short or gave an NA ES), \`n_significant\`,
  \`pct_significant\` (over \`n_valid\`), \`median_es\`, \`iqr_es\`,
  \`sign_consistent\` (proportion of significant cells with ES \> 0).
  For Granger it summarizes x -\> y.

- \`\$robustness_yx\`:

  Granger only: the same fields as \`\$robustness\`, computed from
  \`es_yx\` and \`p_yx\` (y -\> x).

## Details

\*\*Grid construction.\*\* Each vector argument (\`window_sec\`,
\`lag_sec\`, \`increment_pct\`, \`statistic\`, \`surrogate_method\`) is
crossed into a full parameter grid. Seconds are converted to samples per
cell; \`lag_max\` is hard-capped at \`floor(window_size / 2)\` to
preserve statistical reliability. Cells where the series is too short
are silently skipped and appear as \`NA\` rows in the output grid.

\*\*Surrogate reuse.\*\* One surrogate matrix is generated per unique
\`surrogate_method\` and reused across every cell sharing that method;
surrogate cost does not multiply by grid size.

\*\*Effect size polarity.\*\* For WCC and Granger, higher values
indicate stronger synchrony (upper-tail test): \`ES = (obs - null_mean)
/ null_sd\`. For WDTW, lower distance is better (lower-tail test): \`ES
= (null_mean - obs) / null_sd\`. Both have the convention ES \> 0 =
evidence for synchrony.

\*\*Granger direction.\*\* When \`estimator = "wgranger"\`, two sets of
statistics are returned: primary (\`observed\`, \`null_mean\`,
\`null_sd\`, \`es\`, \`p\`) refer to x -\> y; additional columns
\`observed_yx\`, \`null_mean_yx\`, \`null_sd_yx\`, \`es_yx\`, and
\`p_yx\` give y -\> x. \`\$robustness\` summarizes x -\> y and
\`\$robustness_yx\` summarizes y -\> x. The \`print()\`, \`summary()\`,
\`glance()\`, and \`plot()\` methods report x -\> y by default; pass
\`direction = "yx"\` to report y -\> x.

## See also

\[autotune_wcc()\], \[suggest_wcc_params()\],
\[plot.bsync_multiverse()\], \[tidy.bsync_multiverse()\],
\[glance.bsync_multiverse()\];
[`generate_surrogate_iaaft`](https://jmgirard.github.io/bsync/reference/generate_surrogate_iaaft.md)
for the `"iaaft"` method

## Examples

``` r
# \donttest{
# Sweep a seconds-specified window/lag grid and test each cell vs. a null.
# A small surrogate count keeps the example fast; use >= 1000 for reporting.
mv <- synchrony_multiverse(
  x = sim_dyad$x_A,
  y = sim_dyad$x_B,
  estimator = "wcc",
  sample_rate = 80,
  window_sec = c(1, 2, 4),
  lag_sec = 1,
  n_surrogates = 50
)
mv
#> 
#> ── Synchrony Multiverse Analysis (wcc) ─────────────────────────────────────────
#> Specifications: 3 (3 computable)
#> Surrogates per cell: 50
#> Significant (p <= .05): 0 of 3 (0%)
#> Median ES: 1.193 [IQR: 0.491]
#> Sign-consistent (sig. cells): NA%
glance(mv)
#> # A tibble: 1 × 9
#>   estimator n_cells n_valid n_significant pct_significant median_es iqr_es
#>   <chr>       <int>   <int>         <int>           <dbl>     <dbl>  <dbl>
#> 1 wcc             3       3             0               0      1.19  0.491
#> # ℹ 2 more variables: sign_consistent <dbl>, n_surrogates <dbl>
# }
```

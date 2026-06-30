# Select the Best Specification from a Multi-Dyad Multiverse

Applies the gated stability-penalized selection rule to a list of
\`bsync_multiverse\` objects (one per dyad) and returns the winning
cell.

## Usage

``` r
select_specification(mv_list, sig_pct = 0.5, iqr_penalty = 0.5)
```

## Arguments

- mv_list:

  A list of \`bsync_multiverse\` objects, all run with the same
  parameter grid (i.e., all produced by \[synchrony_multiverse()\] with
  identical \`window_sec\`, \`lag_sec\`, \`increment_pct\`,
  \`statistic\`, and \`surrogate_method\` arguments).

- sig_pct:

  Minimum proportion of dyads in which a cell must be significant (p \<
  .05) to pass the detectability gate. Default \`0.5\`.

- iqr_penalty:

  Penalty weight on cross-dyad IQR of ES in the score \`median(ES) -
  iqr_penalty \* IQR(ES)\`. Default \`0.5\`.

## Value

A list with \`best_row\` (one-row tibble from the grid), \`sig_rate\`,
\`median_es\`, \`iqr_es\`, \`score\`, and \`n_gated\` for the selected
cell.

## See also

\[autotune_wcc()\], \[synchrony_multiverse()\]

## Examples

``` r
# \donttest{
# Build one multiverse per dyad, then pick the most robust specification.
# Small surrogate count for a fast example; use >= 1000 for reporting.
mv_list <- lapply(seq_len(3), function(i) {
  synchrony_multiverse(
    x = sim_dyad$x_A,
    y = sim_dyad$x_B,
    estimator = "wcc",
    sample_rate = 80,
    window_sec = c(1, 2, 4),
    lag_sec = 1,
    n_surrogates = 30
  )
})
select_specification(mv_list)
#> Warning: No specification passed the detectability gate (0.5 significance rate). Falling
#> back to highest median ES.
#> $best_row
#> # A tibble: 1 × 15
#>   estimator window_sec lag_sec increment_pct surrogate_method statistic 
#>   <chr>          <dbl>   <dbl>         <dbl> <chr>            <chr>     
#> 1 wcc                2       1           0.1 phase            mean_abs_z
#> # ℹ 9 more variables: window_size <dbl>, lag_max <int>, window_increment <dbl>,
#> #   n_windows <dbl>, observed <dbl>, null_mean <dbl>, null_sd <dbl>, es <dbl>,
#> #   p <dbl>
#> 
#> $sig_rate
#> [1] 0
#> 
#> $median_es
#> [1] 1.329649
#> 
#> $iqr_es
#> [1] 0.324802
#> 
#> $score
#> [1] 1.167248
#> 
#> $n_gated
#> [1] 3
#> 
# }
```

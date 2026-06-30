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

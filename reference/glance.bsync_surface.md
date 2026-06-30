# One-row summary of a bsync_surface object

Returns a single-row tibble with the aggregate statistic(s) and key
settings so multiple estimator runs can be compared with
\[dplyr::bind_rows()\].

## Usage

``` r
# S3 method for class 'bsync_surface'
glance(x, ...)
```

## Arguments

- x:

  A \`bsync_surface\` object (\`wcc_res\`, \`wdtw_res\`, or
  \`wgranger_res\`).

- ...:

  Additional arguments (not used).

## Value

A one-row \[tibble::tibble()\].

## Details

For WCC/WDTW the aggregate is a single named column (\`mean_abs_z\`,
\`peak\`, or \`mean_distance\`). For Granger, \`f_xy\` and \`f_yx\` are
separate columns.

## See also

\[tidy.bsync_surface()\], \[as_tibble.bsync_surface()\]

## Examples

``` r
wcc_res <- wcc(sim_dyad$x_A, sim_dyad$x_B, window_size = 96, lag_max = 10)
# One-row summary: aggregate statistic(s) + key settings
glance(wcc_res)
#> # A tibble: 1 × 7
#>   mean_abs_z n_windows window_size window_increment lag_max lag_increment
#>        <dbl>     <int>       <dbl>            <dbl>   <dbl>         <dbl>
#> 1     0.0900      2285          96                1      10             1
#> # ℹ 1 more variable: statistic <chr>
```

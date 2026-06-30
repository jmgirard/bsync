# Tidy a bsync_surface object into a tibble of per-cell results

Returns one row per cell in \`results_df\`: window position x lag (or
just window position for Granger). Column names match the underlying
estimator.

## Usage

``` r
# S3 method for class 'bsync_surface'
tidy(x, ...)
```

## Arguments

- x:

  A \`bsync_surface\` object (\`wcc_res\`, \`wdtw_res\`, or
  \`wgranger_res\`).

- ...:

  Additional arguments (not used).

## Value

A \[tibble::tibble()\].

## See also

\[glance.bsync_surface()\], \[as_tibble.bsync_surface()\]

## Examples

``` r
wcc_res <- wcc(sim_dyad$x_A, sim_dyad$x_B, window_size = 96, lag_max = 10)
# One row per surface cell (window x lag)
tidy(wcc_res)
#> # A tibble: 47,985 × 3
#>        i   tau     wcc
#>    <dbl> <int>   <dbl>
#>  1    11   -10 -0.116 
#>  2    12   -10 -0.0938
#>  3    13   -10 -0.122 
#>  4    14   -10 -0.111 
#>  5    15   -10 -0.0914
#>  6    16   -10 -0.115 
#>  7    17   -10 -0.131 
#>  8    18   -10 -0.132 
#>  9    19   -10 -0.129 
#> 10    20   -10 -0.126 
#> # ℹ 47,975 more rows
```

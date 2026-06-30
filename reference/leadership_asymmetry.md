# Calculate Leadership Asymmetry Index

Computes a rolling index of leader-follower asymmetry from optimally
picked lags. A value of 1 indicates 'x' leads entirely, -1 indicates 'y'
leads entirely, and 0 indicates equal leading or simultaneous behavior.

## Usage

``` r
leadership_asymmetry(optima_obj, epoch_size = 10, min_valid = 3)
```

## Arguments

- optima_obj:

  An object of class "wcc_optima" or "wdtw_optima".

- epoch_size:

  A positive integer specifying the total width of the centered sliding
  window (in number of optima) used to compute each local asymmetry
  ratio. (default = \`10\`)

- min_valid:

  A positive integer specifying the minimum number of valid (non-\`NA\`)
  optima required inside an epoch to compute the index. Positions with
  fewer valid optima receive \`NA\`. Must be at least 1. (default =
  \`3\`)

## Value

A data frame containing the rolling asymmetry index.

## Details

The index is computed using a \*\*centered sliding window\*\* of
\`epoch_size\` windows. For each position \`t\`, the window spans \`t -
floor(epoch_size / 2)\` to \`t + floor(epoch_size / 2)\` (clamped to the
series boundaries). Within that neighborhood the fraction of positive
vs. negative optimum lags determines the asymmetry score. Positions
where fewer than \`min_valid\` non-\`NA\` optima fall inside the window
receive \`NA\` in the output.

## Examples

``` r
wcc_res <- wcc(sim_dyad$x_A, sim_dyad$x_B, window_size = 96, lag_max = 10)
optima <- pick_optima(wcc_res, L_size = 9)
lai <- leadership_asymmetry(optima, epoch_size = 10)
head(lai)
#> 
#> ── Leadership Asymmetry Index (LAI) ────────────────────────────────────────────
#> Valid Epochs Computed: 6 out of 6
#> Overall Mean LAI: -0.754
#> Overall Dynamics: y predominantly leads
#> Showing the first 5 results:
#>   i x_leads_n y_leads_n asymmetry_index
#>  11         1         5      -0.6666667
#>  12         1         6      -0.7142857
#>  13         1         7      -0.7500000
#>  14         1         8      -0.7777778
#>  15         1         9      -0.8000000
```

# Windowed Granger Causality

Conduct a rolling windowed Granger Causality analysis to determine
dynamic leader-follower relationships between two continuous time
series.

## Usage

``` r
wgranger(x, y, time = NULL, window_size, ar_order = 1, window_increment = 1)
```

## Arguments

- x:

  A numeric vector containing a time series (same length as \`y\`).

- y:

  A numeric vector containing a time series (same length as \`x\`).

- time:

  An optional numeric vector representing the timestamps for the data.

- window_size:

  A positive integer indicating the size of each window.

- ar_order:

  A positive integer specifying the Autoregressive (AR) order. This
  represents the maximum number of lags included in the prediction
  model.

- window_increment:

  A positive integer indicating the window shift increment. (default =
  \`1\`)

## Value

A list object of class "wgranger_res".

## Examples

``` r
# Windowed Granger causality (no lag axis; directional F-statistics)
wgc_res <- wgranger(
  x = sim_dyad$x_A,
  y = sim_dyad$x_B,
  window_size = 96
)
wgc_res
#> 
#> ── Windowed Granger Causality Analysis ─────────────────────────────────────────
#> Total Windows: 2305
#> Window Size: 96
#> AR Order (Lags): 1
```

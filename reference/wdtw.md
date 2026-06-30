# Windowed Dynamic Time Warping

Conduct a windowed dynamic time warping (WDTW) analysis to find the
optimal alignment distance between sliding windows of two time series.

## Usage

``` r
wdtw(
  x,
  y,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1,
  scale_method = c("global", "local", "none"),
  distance_metric = c("L2", "L1")
)
```

## Arguments

- x:

  A numeric vector containing a time series (same length as \`y\`).

- y:

  A numeric vector containing a time series (same length as \`x\`).

- time:

  An optional numeric vector representing the timestamps for the data.
  Must be the same length as \`x\` and \`y\`. Default is \`NULL\`.

- window_size:

  A positive integer indicating the size of each window.

- lag_max:

  A positive integer indicating the maximum lag to try.

- window_increment:

  A positive integer indicating the window shift increment. (default =
  \`1\`)

- lag_increment:

  A positive integer indicating the lag shift increment. (default =
  \`1\`)

- scale_method:

  Character string specifying how to standardize the data. "global"
  standardizes the entire time series before analysis. "local"
  standardizes within each sliding window. "none" applies no scaling.
  (default = \`"global"\`)

- distance_metric:

  Character string specifying the local cost function. "L1" uses
  absolute difference (Manhattan). "L2" uses squared difference
  (Euclidean). (default = \`"L2"\`)

## Value

A list object of class "wdtw_res".

## Examples

``` r
# \donttest{
# Windowed dynamic time warping. DTW is O(window^2) per cell, so this
# example runs on a short subset; use the full series in real analyses.
wdtw_res <- wdtw(
  x = sim_dyad$x_A[1:600],
  y = sim_dyad$x_B[1:600],
  window_size = 96,
  lag_max = 10
)
wdtw_res
#> 
#> ── Windowed Dynamic Time Warping Analysis ──────────────────────────────────────
#> Total Windows: 485
#> Total Lags Tested: 21
#> Window Size: 96
#> Max Lag: 10
#> Scale Method: global
#> Distance Metric: L2
#> Overall Mean Distance: 48.1416
# }
```

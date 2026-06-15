# Calculate Surrogate Windowed Dynamic Time Warping (WDTW)

Calculate Surrogate Windowed Dynamic Time Warping (WDTW)

## Usage

``` r
wdtw_surrogate(
  x,
  y,
  y_surrogates,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1,
  scale_method = c("global", "local", "none"),
  distance_metric = c("L2", "L1"),
  fast_method = FALSE
)
```

## Arguments

- x:

  A numeric vector containing the reference time series.

- y:

  A numeric vector containing the query time series.

- y_surrogates:

  A matrix of surrogate time series for \`y\` (columns are surrogates).

- time:

  An optional numeric vector representing the timestamps for the data.
  Default is \`NULL\`.

- window_size:

  A positive integer indicating the size of the rolling window.

- lag_max:

  A positive integer indicating the maximum lag to try.

- window_increment:

  A positive integer indicating the step size for the rolling window.
  Default is 1.

- lag_increment:

  A positive integer indicating the lag shift increment. Default is 1.

- scale_method:

  Character string specifying how to standardize the data. Default is
  \`"global"\`.

- distance_metric:

  Character string specifying the local cost function. Default is
  \`"L2"\`.

- fast_method:

  Logical. If \`TRUE\`, severely reduces computation time by only
  evaluating the surrogate alignments at a lag of 0. Useful for
  exploratory analysis. Default is \`FALSE\`.

## Value

A list object of class "wdtw_surr".

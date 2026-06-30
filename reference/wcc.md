# Windowed Cross-Correlation

Conduct a windowed cross-correlation analysis

## Usage

``` r
wcc(
  x,
  y,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1,
  na.rm = TRUE,
  statistic = c("mean_abs_z", "peak")
)
```

## Arguments

- x:

  A numeric vector containing a time series (same length as \`y\`).

- y:

  A numeric vector containing a time series (same length as \`x\`).

- time:

  An optional numeric vector representing the timestamps for the data.
  Must be the same length as \`x\` and \`y\`. If provided, the rolling
  window indices will be mapped directly to these timestamps in the
  results, which is highly recommended to maintain accurate timelines if
  edge artifacts were trimmed prior to analysis. Default is \`NULL\`.

- window_size:

  A positive integer indicating the size of each window, i.e., the
  number of elements in each window vector. Boker et al. suggest setting
  the window small enough so that the assumption can be made of little
  change in lead-lag relationships within the number of samples in the
  window but not so small that the reliability for the correlation
  estimate for each sample will be reduced.

- lag_max:

  A positive integer indicating the maximum lag to try between \`x\` and
  \`y\` windows. Boker et al. recommend selecting the greatest interval
  of time separating a behavior from participant \`x\` and a behavior
  from participant \`y\` that would be considered to be of interest.

- window_increment:

  A positive integer indicating the number of samples between successive
  changes in the window for the \`x\` vector. Can be made larger than 1
  to reduce the number of rows in the output matrix. Boker et al.
  recommend setting the window increment as long as possible, but not so
  long that the relation between successive rows in the results matrix
  is lost. (default = \`1\`)

- lag_increment:

  A positive integer indicating the number of samples between successive
  changes in the window for the \`y\` vector (and thus also the interval
  of time separating successive columns in the results matrix). Boker et
  al. recommend setting the lag increment to the longest lag increment
  that still results in related change between successive columns.
  (default = \`1\`)

- na.rm:

  A logical indicating whether to remove missing values from the windows
  when calculating windowed cross-correlations. (default = \`TRUE\`)

- statistic:

  A character string specifying how to aggregate the WCC surface into a
  single number. \`"mean_abs_z"\` (default) takes the mean of absolute
  Fisher's Z values over \*\*all\*\* windows and lags – the SUSY \*mean
  absolute Z\* (Tschacher & Meier, 2020). \`"peak"\` takes the maximum
  absolute Fisher's Z across lags \*\*within each window\*\*, then
  averages those per-window peaks – the rMEA \*best-lag\* convention
  (Boker et al., 2002). Both are larger-is-more-synchrony quantities.
  Pass the same value to \`wcc_surrogate()\` so the null distribution
  matches (see Invariant 2).

## Value

A list object of class "wcc_res" containing the results matrix and
useful summaries of it.

## Examples

``` r
# Windowed cross-correlation on the bundled simulated dyad
wcc_res <- wcc(
  x = sim_dyad$x_A,
  y = sim_dyad$x_B,
  window_size = 96,
  lag_max = 10
)
wcc_res
#> 
#> ── Windowed Cross-Correlation Analysis ─────────────────────────────────────────
#> Total Windows: 2285
#> Total Lags Tested: 21
#> Window Size: 96
#> Max Lag: 10
#> Mean Abs. Fisher's Z: 0.09

# Use the rMEA best-lag "peak" aggregate instead of the SUSY mean
wcc(sim_dyad$x_A, sim_dyad$x_B, window_size = 96, lag_max = 10, statistic = "peak")
#> 
#> ── Windowed Cross-Correlation Analysis ─────────────────────────────────────────
#> Total Windows: 2285
#> Total Lags Tested: 21
#> Window Size: 96
#> Max Lag: 10
#> Mean Peak Abs. Fisher's Z: 0.2369
```

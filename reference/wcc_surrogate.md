# Calculate Surrogate Windowed Cross-Correlations

Calculate Surrogate Windowed Cross-Correlations

## Usage

``` r
wcc_surrogate(
  x,
  y,
  y_surrogates,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1,
  na.rm = TRUE
)
```

## Arguments

- x:

  A numeric vector containing a time series.

- y:

  A numeric vector containing a time series.

- y_surrogates:

  A matrix of surrogate time series for \`y\` (columns are surrogates).

- time:

  An optional numeric vector representing the timestamps for the data.
  Default is \`NULL\`.

- window_size:

  A positive integer indicating the size of each window.

- lag_max:

  A positive integer indicating the maximum lag to try.

- window_increment:

  A positive integer indicating the window shift increment. Default is
  1.

- lag_increment:

  A positive integer indicating the lag shift increment. Default is 1.

- na.rm:

  A logical indicating whether to remove missing values. Default is
  \`TRUE\`.

## Value

A list object of class "wcc_surr".

## Details

The p-value is the proportion of surrogates whose aggregate statistic is
\*\*at least as large as\*\* the observed statistic. The aggregate is
\`mean(abs(Fisher's Z))\` — the same quantity stored in
\`wcc_res\$fisher_z\` — computed identically on both the observed data
and every surrogate, so the null distribution and the observed value are
directly comparable.

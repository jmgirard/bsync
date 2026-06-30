# Calculate Surrogate Windowed Granger Causality

Calculate Surrogate Windowed Granger Causality

## Usage

``` r
wgranger_surrogate(
  x,
  y,
  y_surrogates,
  time = NULL,
  window_size,
  ar_order = 1,
  window_increment = 1
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

  A positive integer indicating the size of the rolling window.

- ar_order:

  A positive integer specifying the Autoregressive (AR) order. Default
  is 1.

- window_increment:

  A positive integer indicating the step size for the rolling window.
  Default is 1.

## Value

A list object of class "wgranger_surr".

## Details

Two p-values are returned: one for x → y and one for y → x. Each is the
proportion of surrogates whose mean F-statistic across windows is \*\*at
least as large as\*\* the corresponding observed mean F-statistic. Both
null distributions are built with the same aggregate (\`mean(f_xy)\` and
\`mean(f_yx)\`) as the observed statistics stored in
\`wgranger_res\$results_df\`, so the null and observed values are
directly comparable.

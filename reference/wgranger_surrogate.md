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

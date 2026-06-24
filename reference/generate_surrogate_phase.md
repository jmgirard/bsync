# Generate Phase-Randomized Surrogates (Fourier Transform)

Generate Phase-Randomized Surrogates (Fourier Transform)

## Usage

``` r
generate_surrogate_phase(y, n_surrogates = 100, trim_odd = FALSE)
```

## Arguments

- y:

  A numeric vector containing a time series.

- n_surrogates:

  Integer specifying the number of surrogates. Default is 100.

- trim_odd:

  Logical. If TRUE, drops the final observation if the time series
  length is odd.

## Value

A matrix where each column is a surrogate time series.

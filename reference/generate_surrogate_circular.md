# Generate Circular Shift Surrogates

Generate Circular Shift Surrogates

## Usage

``` r
generate_surrogate_circular(y, n_surrogates = 100, lag_max = NULL)
```

## Arguments

- y:

  A numeric vector containing a time series.

- n_surrogates:

  Integer specifying the number of surrogates. Default is 100.

- lag_max:

  Optional integer. If provided, ensures shifts are large enough to
  break local autocorrelation.

## Value

A matrix where each column is a surrogate time series.

## Examples

``` r
# Build 100 circular-shift surrogates of one partner's signal
surr <- generate_surrogate_circular(sim_dyad$x_B, n_surrogates = 100)
dim(surr)
#> [1] 2400  100
```

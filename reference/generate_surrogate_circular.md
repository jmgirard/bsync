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

## See also

[`generate_surrogate_phase`](https://jmgirard.github.io/bsync/reference/generate_surrogate_phase.md),
[`generate_surrogate_iaaft`](https://jmgirard.github.io/bsync/reference/generate_surrogate_iaaft.md),
and
[`generate_surrogate_segment`](https://jmgirard.github.io/bsync/reference/generate_surrogate_segment.md)
for the other within-dyad nulls;
[`generate_surrogate_pseudo`](https://jmgirard.github.io/bsync/reference/generate_surrogate_pseudo.md)
for the between-dyad null.

## Examples

``` r
# Build 100 circular-shift surrogates of one partner's signal
surr <- generate_surrogate_circular(sim_dyad$x_B, n_surrogates = 100)
dim(surr)
#> [1] 2400  100
```

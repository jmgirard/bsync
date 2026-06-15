# Impute Missing Values in Continuous Time Series with Metadata

Impute Missing Values in Continuous Time Series with Metadata

## Usage

``` r
impute_ts_gaps(x, method = c("linear", "spline"), maxgap = 5)
```

## Arguments

- x:

  A numeric vector representing the time series.

- method:

  Character string specifying the interpolation method ("linear" or
  "spline").

- maxgap:

  Integer specifying the maximum number of consecutive NAs to impute.
  Gaps larger than this will be left as NA.

## Value

A numeric vector with small gaps imputed, containing an
"imputation_summary" attribute.

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

## Examples

``` r
# Linearly impute gaps up to 5 samples wide; leave longer gaps as NA
x <- c(1, 2, NA, NA, 5, 6, NA, 8, 9, NA, NA, NA, NA, NA, NA, 16)
impute_ts_gaps(x, method = "linear", maxgap = 5)
#> Warning: Found 1 gap exceeding `maxgap` (5); left as NA.
#>  [1]  1  2  3  4  5  6  7  8  9 NA NA NA NA NA NA 16
#> attr(,"imputation_summary")
#> attr(,"imputation_summary")$method
#> [1] "linear"
#> 
#> attr(,"imputation_summary")$maxgap_used
#> [1] 5
#> 
#> attr(,"imputation_summary")$values_imputed
#> [1] 3
#> 
#> attr(,"imputation_summary")$values_left_na
#> [1] 6
#> 
```

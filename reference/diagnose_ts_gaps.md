# Diagnose Missing Data Gaps in a Time Series

Diagnose Missing Data Gaps in a Time Series

## Usage

``` r
diagnose_ts_gaps(x)
```

## Arguments

- x:

  A numeric vector.

## Value

A data frame with summary statistics about missing values.

## Examples

``` r
# Summarize the missing-value runs in a signal with gaps
x <- c(1, 2, NA, NA, 5, 6, NA, 8, 9, NA, NA, NA, 13)
diagnose_ts_gaps(x)
#>   total_obs total_na percent_na total_gaps max_gap_length
#> 1        13        6      46.15          3              3
```

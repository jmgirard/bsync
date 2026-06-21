# Downsample a Time Series Signal via Rolling Aggregation

Reduces the sampling rate of a continuous time series by applying a
rolling aggregation function (median or mean) across non-overlapping
windows.

## Usage

``` r
downsample_signal(x, factor, method = c("median", "mean"), na.rm = TRUE)
```

## Arguments

- x:

  A numeric vector representing the continuous time series signal.

- factor:

  A single positive integer indicating the downsampling factor. For
  example, a factor of 6 reduces 30Hz data to 5Hz.

- method:

  A character string specifying the aggregation method: "median" or
  "mean". Default is "median", which is highly robust to single-frame
  tracking glitches.

- na.rm:

  A logical indicating whether to remove missing values during
  aggregation. Default is \`TRUE\`.

## Value

A numeric vector representing the downsampled time series.

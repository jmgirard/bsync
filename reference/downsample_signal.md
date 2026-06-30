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

## Details

\*\*When to use this function versus \`aggregate_by_time()\`:\*\*

\* Use \`downsample_signal()\` when working with a single, continuous
numeric vector that has guaranteed regular intervals and no missing
frames. That function relies on matrix reshaping and vector math, making
it exceptionally fast for clean, pre-processed data. \* Use
\`aggregate_by_time()\` when working with raw behavioral tracking data
(e.g., OpenFace output) that may contain irregular timestamps, dropped
frames, or missing rows. By binning based on the actual time variable,
this function preserves the true chronological structure of the data and
correctly leaves gaps where tracking was lost. It is also ideal for
processing multiple numeric columns simultaneously.

## Examples

``` r
# Downsample by a factor of 4 (e.g. 80 Hz -> 20 Hz) via the median
ds <- downsample_signal(sim_dyad$x_A, factor = 4)
length(ds)
#> [1] 600
```

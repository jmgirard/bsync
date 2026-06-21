# Aggregate Time Series Data by Time Bins

Efficiently downsamples time series data by aggregating values within
specified time bins. This is a high-level, data frame-based pipeline
function.

## Usage

``` r
aggregate_by_time(
  data,
  time_var,
  bin_width,
  method = c("median", "mean"),
  na.rm = TRUE
)
```

## Arguments

- data:

  A data frame containing the time series data.

- time_var:

  The unquoted name of the column containing time values.

- bin_width:

  A numeric value specifying the width of the time bins. This should be
  in the same units as your time variable (e.g., 0.1 for 100ms bins).

- method:

  A character string specifying the aggregation method: "median" or
  "mean". Default is "median", which is highly robust to single-frame
  tracking glitches.

- na.rm:

  A logical indicating whether to remove missing values when calculating
  the aggregate. Default is \`TRUE\`.

## Value

A new data frame with the downsampled time series. The time variable is
updated to represent the center of each bin, and all non-numeric columns
are dropped.

## Details

\*\*When to use this function versus \`downsample_signal()\`:\*\*

\* Use \`aggregate_by_time()\` when working with raw behavioral tracking
data (e.g., OpenFace output) that may contain irregular timestamps,
dropped frames, or missing rows. By binning based on the actual time
variable, this function preserves the true chronological structure of
the data and correctly leaves gaps where tracking was lost. It is also
ideal for processing multiple numeric columns simultaneously. \* Use
\`downsample_signal()\` when working with a single, continuous numeric
vector that has guaranteed regular intervals and no missing frames. That
function relies on matrix reshaping and vector math, making it
exceptionally fast for clean, pre-processed data.

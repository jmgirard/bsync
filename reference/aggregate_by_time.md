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

## Examples

``` r
# Aggregate the dyad into 0.5-second bins by the median
binned <- aggregate_by_time(sim_dyad, time_var = time, bin_width = 0.5)
head(binned)
#>   time           x_A           y_A         z_A           x_B           y_B
#> 1 0.25 -3.654198e-05  1.258303e-04  0.01755618 -3.369233e-04 -6.753525e-04
#> 2 0.75 -4.049888e-04  6.429244e-05  0.05283492 -9.101645e-05  1.007204e-04
#> 3 1.25 -6.293527e-04  4.994968e-04 -0.09101961  2.880562e-04 -9.402914e-05
#> 4 1.75  7.104979e-04  7.502840e-04 -0.12817631 -2.970780e-04  2.097257e-04
#> 5 2.25  1.100113e-04 -2.720645e-05  0.16352097 -1.801977e-04 -3.462351e-04
#> 6 2.75 -5.013178e-05  2.892379e-04  0.20258062  6.315624e-05  7.454427e-05
#>           z_B
#> 1 -0.01160438
#> 2  0.05875119
#> 3  0.08180088
#> 4 -0.14807526
#> 5 -0.12357705
#> 6  0.24560688
```

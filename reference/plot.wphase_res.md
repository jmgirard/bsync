# Plot wphase_res object

Plot wphase_res object

## Usage

``` r
# S3 method for class 'wphase_res'
plot(
  x,
  time_step = 1,
  color_low = "#F7F7F7",
  color_high = "#2166AC",
  show_zero_lag = TRUE,
  zero_line_color = "black",
  ...
)
```

## Arguments

- x:

  An object of class "wphase_res".

- time_step:

  A numeric value specifying the duration of each index. If not 1, axes
  will be converted from raw indices to time units. Default is 1.

- color_low:

  Character string specifying the color for a PLV of 0. Default is
  "#F7F7F7" (Off-white).

- color_high:

  Character string specifying the color for a PLV of 1. Default is
  "#2166AC" (Deep Blue).

- show_zero_lag:

  Logical indicating whether to draw a vertical line at lag = 0. Default
  is \`TRUE\`.

- zero_line_color:

  Character string specifying the color of the zero-lag line. Default is
  "black".

- ...:

  Additional arguments (not used).

## Value

A \`ggplot2\` plot object.

## Examples

``` r
res <- wphase(sim_dyad$z_A, sim_dyad$z_B, window_size = 96, lag_max = 10)
plot(res)
```

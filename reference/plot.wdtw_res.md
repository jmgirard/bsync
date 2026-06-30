# Plot wdtw_res object

Plot wdtw_res object

## Usage

``` r
# S3 method for class 'wdtw_res'
plot(x, time_step = 1, show_zero_lag = TRUE, zero_line_color = "black", ...)
```

## Arguments

- x:

  An object of class "wdtw_res".

- time_step:

  A numeric value specifying the duration of each index. If not 1, axes
  will be converted from raw indices to time units. Default is 1.

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

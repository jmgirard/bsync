# Plot Surface with Optima Overlay

Plot Surface with Optima Overlay

## Usage

``` r
plot_optima_overlay(
  surface_obj,
  optima_df,
  time_step = 1,
  line_color = "black",
  point_fill = "black",
  point_stroke = "white",
  show_zero_lag = TRUE,
  zero_line_color = "black",
  ...
)
```

## Arguments

- surface_obj:

  An object of class "wcc_res" or "wdtw_res".

- optima_df:

  A data frame of class "wcc_optima" or "wdtw_optima".

- time_step:

  A numeric value specifying the duration of each index. Default is 1.

- line_color:

  Character string specifying the color of the connecting line. Default
  is "black".

- point_fill:

  Character string specifying the inner fill of the optima points.
  Default is "black".

- point_stroke:

  Character string specifying the outer outline of the optima points.
  Default is "white".

- show_zero_lag:

  Logical indicating whether to draw a vertical line at lag = 0. Default
  is \`TRUE\`.

- zero_line_color:

  Character string specifying the color of the zero-lag line. Default is
  "black".

- ...:

  Additional arguments passed to the underlying plot method (e.g.,
  custom colors for WCC).

## Value

A \`ggplot2\` plot object.

## Examples

``` r
wcc_res <- wcc(sim_dyad$x_A, sim_dyad$x_B, window_size = 96, lag_max = 10)
optima <- pick_optima(wcc_res, L_size = 9)
plot_optima_overlay(wcc_res, optima)
```

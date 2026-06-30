# Calculate 1D Speed

Calculates the absolute speed along a single axis (e.g., x-axis) using
finite difference methods.

## Usage

``` r
calc_speed_1d(
  t,
  x,
  n = 1,
  method = c("central", "forward", "backward"),
  fill_edges = TRUE
)
```

## Arguments

- t:

  A numeric vector representing time.

- x:

  A numeric vector of coordinates.

- n:

  An integer specifying the step size for the difference calculations.
  Default is 1.

- method:

  A character string specifying the method to use: "central", "forward",
  or "backward". Default is "central".

- fill_edges:

  A logical indicating whether to automatically use forward and backward
  differences to estimate speeds at the boundaries when \`method =
  "central"\`. Default is \`TRUE\`.

## Value

A numeric vector of speeds the same length as the input vectors.

## Examples

``` r
# Unsigned 1D speed (absolute rate of change)
s <- calc_speed_1d(t = sim_dyad$time, x = sim_dyad$x_A)
head(s)
#> [1] 0.255938088 0.030495365 0.079562193 0.064443356 0.194426175 0.005478291
```

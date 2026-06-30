# Calculate 1D Velocity

Calculates the velocity (rate of change including direction) along a
single axis using finite difference methods.

## Usage

``` r
calc_velocity_1d(
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
  differences to estimate velocities at the boundaries when \`method =
  "central"\`. Default is \`TRUE\`.

## Value

A numeric vector of velocities the same length as the input vectors.

## Examples

``` r
# Signed 1D velocity of one coordinate over time
v <- calc_velocity_1d(t = sim_dyad$time, x = sim_dyad$x_A)
head(v)
#> [1] -0.255938088 -0.030495365  0.079562193 -0.064443356 -0.194426175
#> [6] -0.005478291
```

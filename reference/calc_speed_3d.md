# Calculate 3D Speed

Calculates the magnitude of the 3D velocity vector (speed) using finite
difference methods.

## Usage

``` r
calc_speed_3d(
  t,
  x,
  y,
  z,
  n = 1,
  method = c("central", "forward", "backward"),
  fill_edges = TRUE
)
```

## Arguments

- t:

  A numeric vector representing time.

- x:

  A numeric vector of x-coordinates.

- y:

  A numeric vector of y-coordinates.

- z:

  A numeric vector of z-coordinates.

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
# 3D speed from x/y/z coordinate channels
s <- calc_speed_3d(
  t = sim_dyad$time, x = sim_dyad$x_A, y = sim_dyad$y_A, z = sim_dyad$z_A
)
head(s)
#> [1] 0.40447985 0.04907835 0.25211593 0.18896803 0.27919165 0.26524432
```

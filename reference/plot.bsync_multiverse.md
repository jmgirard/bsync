# Plot a synchrony multiverse specification curve

Draws a Simonsohn-style specification curve for a `bsync_multiverse`
object. The top panel shows effect sizes sorted from smallest to
largest, with significant cells (p \<= .05) highlighted. The bottom
panel is a choice dashboard showing which analytic choices each
specification used.

## Usage

``` r
# S3 method for class 'bsync_multiverse'
plot(
  x,
  sig_color = "#2166AC",
  insig_color = "grey60",
  active_color = "#2166AC",
  point_size = 1.5,
  top_frac = 0.55,
  direction = c("xy", "yx"),
  ...
)
```

## Arguments

- x:

  A `bsync_multiverse` object.

- sig_color:

  Color for significant cells (p \<= .05). Default: `"#2166AC"`.

- insig_color:

  Color for non-significant cells. Default: `"grey60"`.

- active_color:

  Fill for active choice tiles. Default: `"#2166AC"`.

- point_size:

  Size of ES points. Default: `1.5`.

- top_frac:

  Fraction of plot height allocated to the ES panel. Default: `0.55`.

- direction:

  For a Granger result (`estimator = "wgranger"`), the direction to
  plot: `"xy"` (x -\> y, the default) or `"yx"` (y -\> x). The title
  names the direction. Other estimators have one direction only, so
  `"yx"` is an error for them. Pass it by name (`direction = "yx"`): it
  comes after the color and size arguments.

- ...:

  Additional arguments (not used).

## Value

Returns `x` invisibly; draws to the active graphics device.

## See also

[`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md),
[`tidy.bsync_multiverse()`](https://jmgirard.github.io/bsync/reference/tidy.bsync_multiverse.md),
[`glance.bsync_multiverse()`](https://jmgirard.github.io/bsync/reference/glance.bsync_multiverse.md)

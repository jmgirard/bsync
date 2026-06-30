# Plot a synchrony multiverse specification curve

Draws a Simonsohn-style specification curve for a \`bsync_multiverse\`
object. The top panel shows effect sizes sorted from smallest to
largest, with significant cells (p \< .05) highlighted. The bottom panel
is a choice dashboard showing which analytic choices each specification
used.

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
  ...
)
```

## Arguments

- x:

  A \`bsync_multiverse\` object.

- sig_color:

  Color for significant cells (p \< .05). Default: \`"#2166AC"\`.

- insig_color:

  Color for non-significant cells. Default: \`"grey60"\`.

- active_color:

  Fill for active choice tiles. Default: \`"#2166AC"\`.

- point_size:

  Size of ES points. Default: \`1.5\`.

- top_frac:

  Fraction of plot height allocated to the ES panel. Default: \`0.55\`.

- ...:

  Additional arguments (not used).

## Value

Returns \`x\` invisibly; draws to the active graphics device.

## See also

\[synchrony_multiverse()\], \[tidy.bsync_multiverse()\],
\[glance.bsync_multiverse()\]

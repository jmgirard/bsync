# One-row robustness summary of a bsync_multiverse object

Returns a single-row tibble summarising robustness across the
specification curve: the number of cells, significance rate, median
effect size, IQR, and sign-consistency.

## Usage

``` r
# S3 method for class 'bsync_multiverse'
glance(x, ...)
```

## Arguments

- x:

  A \`bsync_multiverse\` object from \[synchrony_multiverse()\].

- ...:

  Additional arguments (not used).

## Value

A one-row \[tibble::tibble()\].

## See also

\[tidy.bsync_multiverse()\], \[as_tibble.bsync_multiverse()\],
\[synchrony_multiverse()\]

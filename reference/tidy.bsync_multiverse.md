# Tidy a bsync_multiverse object into the specification grid

Returns the parameter grid as a tibble: one row per specification cell
with all analytic choices and the resulting effect size, p-value, and
null statistics.

## Usage

``` r
# S3 method for class 'bsync_multiverse'
tidy(x, ...)
```

## Arguments

- x:

  A \`bsync_multiverse\` object from \[synchrony_multiverse()\].

- ...:

  Additional arguments (not used).

## Value

A \[tibble::tibble()\] (the \`\$grid\` slot).

## See also

\[glance.bsync_multiverse()\], \[as_tibble.bsync_multiverse()\],
\[synchrony_multiverse()\]

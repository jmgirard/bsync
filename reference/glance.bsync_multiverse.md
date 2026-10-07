# One-row robustness summary of a bsync_multiverse object

Returns a single-row tibble summarising robustness across the
specification curve: the number of cells, significance rate, median
effect size, IQR, and sign-consistency.

## Usage

``` r
# S3 method for class 'bsync_multiverse'
glance(x, direction = c("xy", "yx"), ...)
```

## Arguments

- x:

  A `bsync_multiverse` object from
  [`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md).

- direction:

  For a Granger result (`estimator = "wgranger"`), the direction to
  summarize: `"xy"` (x -\> y, the default) or `"yx"` (y -\> x). Other
  estimators have one direction only, so `"yx"` is an error for them.

- ...:

  Additional arguments (not used).

## Value

A one-row
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html).
For a Granger result, a `direction` column after `estimator` names the
direction summarized.

## See also

[`tidy.bsync_multiverse()`](https://jmgirard.github.io/bsync/reference/tidy.bsync_multiverse.md),
[`as_tibble.bsync_multiverse()`](https://jmgirard.github.io/bsync/reference/as_tibble.bsync_multiverse.md),
[`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md)

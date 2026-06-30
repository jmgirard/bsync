# One-row summary of a bsync_surface object

Returns a single-row tibble with the aggregate statistic(s) and key
settings so multiple estimator runs can be compared with
\[dplyr::bind_rows()\].

## Usage

``` r
# S3 method for class 'bsync_surface'
glance(x, ...)
```

## Arguments

- x:

  A \`bsync_surface\` object (\`wcc_res\`, \`wdtw_res\`, or
  \`wgranger_res\`).

- ...:

  Additional arguments (not used).

## Value

A one-row \[tibble::tibble()\].

## Details

For WCC/WDTW the aggregate is a single named column (\`mean_abs_z\`,
\`peak\`, or \`mean_distance\`). For Granger, \`f_xy\` and \`f_yx\` are
separate columns.

## See also

\[tidy.bsync_surface()\], \[as_tibble.bsync_surface()\]

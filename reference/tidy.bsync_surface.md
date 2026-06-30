# Tidy a bsync_surface object into a tibble of per-cell results

Returns one row per cell in \`results_df\`: window position x lag (or
just window position for Granger). Column names match the underlying
estimator.

## Usage

``` r
# S3 method for class 'bsync_surface'
tidy(x, ...)
```

## Arguments

- x:

  A \`bsync_surface\` object (\`wcc_res\`, \`wdtw_res\`, or
  \`wgranger_res\`).

- ...:

  Additional arguments (not used).

## Value

A \[tibble::tibble()\].

## See also

\[glance.bsync_surface()\], \[as_tibble.bsync_surface()\]

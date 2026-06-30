# Plot wgranger_res object

Plot wgranger_res object

## Usage

``` r
# S3 method for class 'wgranger_res'
plot(x, time_step = 1, metric = c("F", "p"), smooth = FALSE, ...)
```

## Arguments

- x:

  An object of class "wgranger_res".

- time_step:

  A numeric value specifying the duration of each index. Default is 1.

- metric:

  Character string specifying which metric to plot. Options are \`"F"\`
  for F-Statistics or \`"p"\` for negative log10 p-values. Default is
  \`"F"\`.

- smooth:

  Logical indicating whether to apply loess smoothing to the lines.
  Default is \`FALSE\`.

- ...:

  Additional arguments (not used).

## Value

A \`ggplot2\` plot object.

# Plot bsync_lai object

Plot bsync_lai object

## Usage

``` r
# S3 method for class 'bsync_lai'
plot(x, time_step = 1, line_color = "#2166AC", smooth = FALSE, ...)
```

## Arguments

- x:

  An object of class "bsync_lai".

- time_step:

  A numeric value specifying the duration of each index. Default is 1.

- line_color:

  Character string specifying the color of the LAI line. Default is
  "black".

- smooth:

  Logical indicating whether to add a loess smoothing line. Default is
  \`FALSE\`.

- ...:

  Additional arguments (not used).

## Value

A \`ggplot2\` plot object.

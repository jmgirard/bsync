# Trim Edge Effects from Data

Removes or masks a specified number of observations from the beginning
and end of a vector or data frame. This is mathematically required after
applying symmetric rolling filters (like Savitzky-Golay) to remove
boundary artifacts.

## Usage

``` r
trim_edges(x, trim_length, pad_na = FALSE)
```

## Arguments

- x:

  A numeric vector, matrix, or data frame.

- trim_length:

  An integer specifying the number of observations to mask from both
  ends. Best practice: For a Savitzky-Golay filter, this must exactly
  equal (window - 1) / 2.

- pad_na:

  A logical indicating whether to replace the trimmed edges with \`NA\`
  instead of dropping them. Set to \`TRUE\` when using inside
  \`dplyr::mutate()\` to preserve the original vector length. Default is
  \`FALSE\`.

## Value

An object of the same class as \`x\` with the edges removed or masked.

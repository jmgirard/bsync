# Smooth a Time Series Signal

Applies a smoothing filter to a numeric vector. Smoothing is highly
recommended prior to calculating velocity or running windowed
cross-correlation (WCC) to reduce high-frequency noise and prevent
spurious correlations.

## Usage

``` r
smooth_signal(
  x,
  method = c("sgolay", "moving_average", "butterworth"),
  window = 5,
  sg_order = 3,
  bw_cutoff = 0.1,
  bw_order = 2,
  lower_bound = NULL,
  upper_bound = NULL
)
```

## Arguments

- x:

  A numeric vector representing the signal to be smoothed.

- method:

  A character string specifying the smoothing method: "moving_average",
  "sgolay" (Savitzky-Golay), or "butterworth". Default is "sgolay".

- window:

  An integer specifying the window size (number of data points) for the
  "moving_average" and "sgolay" methods. Must be an odd number for
  "sgolay". Default is 5.

- sg_order:

  An integer specifying the polynomial order for the Savitzky-Golay
  filter. Must be less than \`window\`. Default is 3.

- bw_cutoff:

  A numeric value between 0 and 1 specifying the normalized cutoff
  frequency for the Butterworth filter. Default is 0.1.

- bw_order:

  An integer specifying the order of the Butterworth filter. Default is
  2.

- lower_bound:

  Numeric. If provided, smoothed values below this are clamped to this
  value.

- upper_bound:

  Numeric. If provided, smoothed values above this are clamped to this
  value.

## Value

A numeric vector containing the smoothed signal, of the same length as
\`x\`.

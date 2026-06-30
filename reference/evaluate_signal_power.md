# Evaluate Signal Power and Suggest Downsampling Rate

Calculates the Power Spectral Density (PSD) of continuous time series
using Welch's method. It determines the frequency below which a
specified proportion of the total signal power is captured and
recommends an optimal integer downsampling factor that yields a clean
sampling rate.

## Usage

``` r
evaluate_signal_power(x, sample_rate, threshold = 0.95, quiet = FALSE)
```

## Arguments

- x:

  A numeric vector, a list of numeric vectors, or a data frame. If a
  data frame is provided, all numeric columns will be evaluated.

- sample_rate:

  A single positive number indicating the sampling rate in Hertz.

- threshold:

  A single numeric value between 0 and 1 indicating the cumulative
  proportion of power to capture. Default is 0.95 (95 percent).

- quiet:

  A logical indicating whether to suppress console output. Default is
  \`FALSE\`.

## Value

A list of class \`"signal_power_res"\` containing the calculated cutoff
frequencies, the recommended integer downsampling factor, and the
resulting target frequency. Call \`plot()\` on the result to visualize
cumulative power.

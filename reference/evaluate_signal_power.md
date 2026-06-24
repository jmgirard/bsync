# Evaluate Signal Power and Suggest Downsampling Rate

Calculates the Power Spectral Density (PSD) of continuous time series
using Welch's method. It determines the frequency below which a
specified proportion of the total signal power is captured and
recommends an optimal integer downsampling factor that yields a clean
sampling rate.

## Usage

``` r
evaluate_signal_power(x, sample_rate, threshold = 0.95, plot = TRUE)
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

- plot:

  A logical indicating whether to return a cumulative power plot.
  Default is \`TRUE\`.

## Value

A list containing the calculated cutoff frequencies, the recommended
integer downsampling factor, the resulting target frequency, and
optionally a \`ggplot\` object.

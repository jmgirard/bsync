# Suggest WCC Hyperparameters

Derives principled starting values for Windowed Cross-Correlation
parameters from the \*\*measured signal\*\* via PSD-based dominant
timescale estimation.

## Usage

``` r
suggest_wcc_params(
  x,
  y,
  sample_rate,
  event_duration_sec = NULL,
  max_delay_sec = 3,
  overlap_pct = 0.5,
  min_window_samples = 20L
)
```

## Arguments

- x:

  Numeric vector; the reference time series.

- y:

  Numeric vector; the query time series (same length as \`x\`).

- sample_rate:

  A numeric value indicating the sampling rate in Hertz.

- event_duration_sec:

  Optional numeric override for the dominant behavioral cycle duration
  in seconds. Default \`NULL\` derives it from the PSD of \`x\` and
  \`y\` via \[evaluate_signal_power()\].

- max_delay_sec:

  The maximum plausible reaction time between participants in seconds,
  used to set the initial \`lag_max\`. Default is \`3\`.

- overlap_pct:

  The desired proportion of overlap between consecutive windows (0-1).
  Default is \`0.5\` (50% overlap).

- min_window_samples:

  Minimum number of samples required in a window for a stable
  correlation. Default is \`20\`.

## Value

A named list with \`window_size\`, \`lag_max\`, \`window_increment\`,
and \`lag_increment\`, ready to pass to \[wcc()\].

## Details

\*\*Dominant timescale.\*\* When \`event_duration_sec\` is \`NULL\`
(default), the function estimates the dominant behavioral cycle from the
measured signal via \[evaluate_signal_power()\]: \`event_duration_sec =
1 / primary_cutoff_freq\`. Pass a numeric value to override with your
own theoretical estimate.

\*\*Window size (4-cycles heuristic).\*\* \`window_size =
round(event_duration_sec \* 4 \* sample_rate)\` (Boker et al., 2002).
Four cycles yields a stable within-window correlation across the range
of lead-lag relationships; two cycles (the Nyquist minimum) is too
noisy.

\*\*Hard constraints applied and reported:\*\*

- \`lag_max \<= floor(window_size / 2)\` – the SUSY reliability
  constraint (\`segment \>= 2\*maxlag\`); beyond this the lagged windows
  share fewer than half their samples.

- \`window_size \<= floor(series_length / 2)\` – ensures at least two
  non-overlapping windows fit in the series.

- \`window_size \>= min_window_samples\` – a minimum-samples floor for a
  stable correlation estimate.

All violations produce informative messages, not silent changes.

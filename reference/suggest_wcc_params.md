# Suggest WCC Hyperparameters

Derives principled starting values for Windowed Cross-Correlation
parameters from the \*\*measured signal\*\* via PSD-based timescale
estimation.

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

  Optional numeric override for the characteristic event timescale in
  seconds. Default \`NULL\` derives it from the power-cutoff frequency
  of \`x\` and \`y\` via \[evaluate_signal_power()\].

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

\*\*Signal timescale (PSD power cutoff).\*\* When \`event_duration_sec\`
is \`NULL\` (default), the function estimates a characteristic timescale
from the measured signal via \[evaluate_signal_power()\]:
\`event_duration_sec = 1 / primary_cutoff_freq\`, where
\`primary_cutoff_freq\` is the frequency below which \`threshold\`
(default 95 frequency – the fastest behaviorally relevant timescale – is
used deliberately rather than the single largest spectral peak: raw
movement spectra are dominated by low-frequency / DC components (a slow
drift can hold most of the power), so a spectral-peak estimate would
collapse toward 0 Hz and an unusably long window. The power cutoff is
robust to that. Pass a numeric value to override with your own
theoretical estimate.

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

## Examples

``` r
# Derive starting parameters from the signal's own power-cutoff timescale.
# (The default max delay exceeds window/2 here, so the SUSY rule caps
# lag_max -- an example of the loud, non-destructive constraint messages.)
params <- suggest_wcc_params(
  x = sim_dyad$z_A,
  y = sim_dyad$z_B,
  sample_rate = 80
)
#> ℹ PSD power-cutoff timescale: 0.8 s (cutoff 1.25 Hz).
#> Warning: Requested `max_delay_sec` (3 s = 240 samples) exceeds window_size/2 (128).
#> ℹ Capping `lag_max` at 128 (= 1.6 s).
#> 
#> ── Suggested WCC Parameters ────────────────────────────────────────────────────
#> window_size: 256 (3.2 s)
#> lag_max: 128 (1.6 s)
#> window_increment: 128 (50% overlap)
#> lag_increment: 1
params
#> $window_size
#> [1] 256
#> 
#> $lag_max
#> [1] 128
#> 
#> $window_increment
#> [1] 128
#> 
#> $lag_increment
#> [1] 1
#> 
```

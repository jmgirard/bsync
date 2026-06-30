# Suggest WCC Hyperparameters

Calculates principled starting values for Windowed Cross-Correlation
parameters based on the sampling rate of the data and the theoretical
timing of the behaviors.

## Usage

``` r
suggest_wcc_params(
  sample_rate,
  event_duration_sec = 2,
  max_delay_sec = 3,
  overlap_pct = 0.5
)
```

## Arguments

- sample_rate:

  A numeric value indicating the sampling rate in Hertz (frames per
  second).

- event_duration_sec:

  The expected duration of a single behavioral event in seconds. Used as
  the basis for the 4-cycles-per-window heuristic: \`window_size =
  round(event_duration_sec \* 4 \* sample_rate)\`. Default is 2 (typical
  for brief conversational gestures).

- max_delay_sec:

  The maximum plausible reaction time between participants in seconds.
  Default is 3.

- overlap_pct:

  The desired percentage of overlap between consecutive time windows.
  Default is 0.5 (50 percent overlap).

## Value

A list of recommended parameters ready to be passed to \`wcc()\`.

## Details

The \`window_size\` is derived from the \*\*4-cycles-per-window
heuristic\*\*: a window should span approximately 4 full cycles of the
behavior of interest so that the within-window correlation estimate is
stable across a range of lead–lag relationships (Boker et al., 2002).
Concretely, \`window_size = round(event_duration_sec \* 4 \*
sample_rate)\`. Four cycles ensures enough oscillation to estimate a
reliable correlation, whereas two cycles (the Nyquist minimum) would
leave the estimate too noisy.

The \`lag_max\` is capped at half the \`window_size\` when the requested
\`max_delay_sec\` would exceed it; beyond that point the lagged window
and the reference window share fewer than half their samples, severely
degrading reliability.

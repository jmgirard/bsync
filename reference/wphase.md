# Windowed Phase Synchrony

Conduct a windowed phase-synchrony analysis. Instantaneous phase is
extracted once per series with the analytic signal
([`gsignal::hilbert()`](https://rdrr.io/pkg/gsignal/man/hilbert.html));
each window-by-lag cell then summarizes the phase differences into a
phase-locking value (PLV) and a mean relative phase.

## Usage

``` r
wphase(
  x,
  y,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1
)
```

## Arguments

- x:

  A numeric vector containing a time series (same length as \`y\`).

- y:

  A numeric vector containing a time series (same length as \`x\`).

- time:

  An optional numeric vector of timestamps for the data, same length as
  \`x\` and \`y\`; if provided, window indices in the results are mapped
  to these timestamps (and it must not contain missing values). Default
  is \`NULL\`.

- window_size:

  A positive integer indicating the size of each window in samples.

- lag_max:

  A positive integer indicating the maximum lag (in samples) to try
  between the \`x\` and \`y\` windows.

- window_increment:

  A positive integer indicating the number of samples between successive
  windows. (default = \`1\`)

- lag_increment:

  A positive integer indicating the number of samples between successive
  lags. (default = \`1\`)

## Value

A list object of class \`"wphase_res"\` (a \`bsync_surface\`) with
\`results_df\` (columns \`i\`, \`tau\`, \`plv\`, \`rel_phase\`),
\`aggregate\` (\`mean_plv\`), and \`settings\`.

## Details

The PLV is the modulus of the mean unit phasor of the phase differences
— the phase-locking value of Lachaux, Rodriguez, Martinerie, and Varela
(1999). Lachaux et al. average across trials at a fixed latency;
\`wphase()\` is a single-trial continuous-data adaptation that averages
across the time samples within each sliding window (at each lag offset).
PLV is 1 when the phase difference is constant across the window and
near 0 when it is uniformly scattered. \`rel_phase\` is the circular
mean of the phase differences (radians); positive values mean \`x\`'s
phase is ahead of \`y\`'s.

\*\*Narrowband caveat (stated loudly, per the source method):\*\*
instantaneous phase is well-defined for narrowband signals. Lachaux et
al. band-pass filter before phase extraction; \`wphase()\` leaves
filtering to you (e.g. \[smooth_signal()\] or an external band-pass)
rather than silently changing the basis. Beware that this failure mode
is \*\*anti-conservative\*\*: broadband or low-frequency-dominated input
\*inflates\* PLV toward 1 (two independent AR(0.98) series measure a
mean PLV near 0.6 at \`window_size = 96\`, against roughly 0.12 for
white noise), because slowly drifting phases stay aligned across a short
window. High PLV on unfiltered signals is not evidence of synchrony;
band-pass first and use \[wphase_surrogate()\] for significance.

\*\*NA policy:\*\* the analytic signal is FFT-based, so a single \`NA\`
corrupts every phase estimate; unlike \[wcc()\], there is no coherent
per-window \`na.rm\`. \`wphase()\` therefore aborts on \`NA\`-containing
input — impute with \[impute_ts_gaps()\] or trim with \[trim_edges()\]
(\`pad_na = FALSE\`) first.

## Examples

``` r
# Windowed phase synchrony on the bundled simulated dyad
wphase_res <- wphase(
  x = sim_dyad$x_A,
  y = sim_dyad$x_B,
  window_size = 96,
  lag_max = 10
)
wphase_res
#> 
#> ── Windowed Phase Synchrony Analysis ───────────────────────────────────────────
#> Total Windows: 2285
#> Total Lags Tested: 21
#> Window Size: 96
#> Max Lag: 10
#> Mean PLV: 0.1196
```

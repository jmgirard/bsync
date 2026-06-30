# Package index

## Estimators

Windowed estimators of nonstationary lead-lag synchrony. All three
return a `bsync_surface` object with a shared
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- [`wcc()`](https://jmgirard.github.io/bsync/reference/wcc.md) :
  Windowed Cross-Correlation
- [`wdtw()`](https://jmgirard.github.io/bsync/reference/wdtw.md) :
  Windowed Dynamic Time Warping
- [`wgranger()`](https://jmgirard.github.io/bsync/reference/wgranger.md)
  : Windowed Granger Causality
- [`print(`*`<wcc_res>`*`)`](https://jmgirard.github.io/bsync/reference/print.wcc_res.md)
  : Print method for wcc_res objects
- [`print(`*`<wdtw_res>`*`)`](https://jmgirard.github.io/bsync/reference/print.wdtw_res.md)
  : Print method for wdtw_res objects
- [`print(`*`<wgranger_res>`*`)`](https://jmgirard.github.io/bsync/reference/print.wgranger_res.md)
  : Print method for wgranger_res objects
- [`summary(`*`<wcc_res>`*`)`](https://jmgirard.github.io/bsync/reference/summary.wcc_res.md)
  : Summary method for wcc_res objects
- [`summary(`*`<wgranger_res>`*`)`](https://jmgirard.github.io/bsync/reference/summary.wgranger_res.md)
  : Summary method for wgranger_res objects
- [`plot(`*`<wcc_res>`*`)`](https://jmgirard.github.io/bsync/reference/plot.wcc_res.md)
  : Plot wcc_res object
- [`plot(`*`<wdtw_res>`*`)`](https://jmgirard.github.io/bsync/reference/plot.wdtw_res.md)
  : Plot wdtw_res object
- [`plot(`*`<wgranger_res>`*`)`](https://jmgirard.github.io/bsync/reference/plot.wgranger_res.md)
  : Plot wgranger_res object

## Surrogate testing

Generate null distributions (pseudo-synchrony) and compare observed
statistics against them. Two-step pipeline: generate a surrogate matrix,
then pass it to the matched surrogate wrapper.

- [`generate_surrogate_circular()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_circular.md)
  : Generate Circular Shift Surrogates
- [`generate_surrogate_phase()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_phase.md)
  : Generate Phase-Randomized Surrogates (Fourier Transform)
- [`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md)
  : Calculate Surrogate Windowed Cross-Correlations
- [`wdtw_surrogate()`](https://jmgirard.github.io/bsync/reference/wdtw_surrogate.md)
  : Calculate Surrogate Windowed Dynamic Time Warping (WDTW)
- [`wgranger_surrogate()`](https://jmgirard.github.io/bsync/reference/wgranger_surrogate.md)
  : Calculate Surrogate Windowed Granger Causality
- [`print(`*`<wcc_surr>`*`)`](https://jmgirard.github.io/bsync/reference/print.wcc_surr.md)
  : Print method for wcc_surr objects
- [`print(`*`<wdtw_surr>`*`)`](https://jmgirard.github.io/bsync/reference/print.wdtw_surr.md)
  : Print method for wdtw_surr objects
- [`print(`*`<wgranger_surr>`*`)`](https://jmgirard.github.io/bsync/reference/print.wgranger_surr.md)
  : Print method for wgranger_surr objects

## Optima extraction & leadership

Extract lag optima from an estimator surface and derive a continuous
Leadership Asymmetry Index (LAI).

- [`pick_optima()`](https://jmgirard.github.io/bsync/reference/pick_optima.md)
  : Find Optimum (Peak or Valley) in Windowed Analyses
- [`leadership_asymmetry()`](https://jmgirard.github.io/bsync/reference/leadership_asymmetry.md)
  : Calculate Leadership Asymmetry Index
- [`plot_optima_overlay()`](https://jmgirard.github.io/bsync/reference/plot_optima_overlay.md)
  : Plot Surface with Optima Overlay
- [`print(`*`<wcc_optima>`*`)`](https://jmgirard.github.io/bsync/reference/print.wcc_optima.md)
  : Print method for wcc_optima objects
- [`print(`*`<wdtw_optima>`*`)`](https://jmgirard.github.io/bsync/reference/print.wdtw_optima.md)
  : Print method for wdtw_optima objects
- [`summary(`*`<wcc_optima>`*`)`](https://jmgirard.github.io/bsync/reference/summary.wcc_optima.md)
  : Summary method for wcc_optima objects
- [`summary(`*`<wdtw_optima>`*`)`](https://jmgirard.github.io/bsync/reference/summary.wdtw_optima.md)
  : Summary method for wdtw_optima objects
- [`print(`*`<bsync_lai>`*`)`](https://jmgirard.github.io/bsync/reference/print.bsync_lai.md)
  : Print method for bsync_lai objects
- [`plot(`*`<bsync_lai>`*`)`](https://jmgirard.github.io/bsync/reference/plot.bsync_lai.md)
  : Plot bsync_lai object

## Parameter guidance

Tools for choosing and evaluating analysis hyperparameters. Start with
[`suggest_wcc_params()`](https://jmgirard.github.io/bsync/reference/suggest_wcc_params.md)
for a single dyad, use
[`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md)
to visualize the specification curve, and
[`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md)
for multi-dyad datasets.
[`select_specification()`](https://jmgirard.github.io/bsync/reference/select_specification.md)
is the underlying selection rule for advanced use.

- [`suggest_wcc_params()`](https://jmgirard.github.io/bsync/reference/suggest_wcc_params.md)
  : Suggest WCC Hyperparameters
- [`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md)
  : Synchrony Multiverse Analysis
- [`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md)
  : Auto-Tune WCC Parameters for a Multi-Dyad Dataset
- [`select_specification()`](https://jmgirard.github.io/bsync/reference/select_specification.md)
  : Select the Best Specification from a Multi-Dyad Multiverse
- [`print(`*`<bsync_autotune>`*`)`](https://jmgirard.github.io/bsync/reference/print.bsync_autotune.md)
  : Print method for bsync_autotune objects
- [`print(`*`<bsync_multiverse>`*`)`](https://jmgirard.github.io/bsync/reference/print.bsync_multiverse.md)
  : Print method for bsync_multiverse objects
- [`summary(`*`<bsync_multiverse>`*`)`](https://jmgirard.github.io/bsync/reference/summary.bsync_multiverse.md)
  : Summary method for bsync_multiverse objects
- [`plot(`*`<bsync_multiverse>`*`)`](https://jmgirard.github.io/bsync/reference/plot.bsync_multiverse.md)
  : Plot a synchrony multiverse specification curve

## Preprocessing — signal & kinematics

Smooth, trim, and compute kinematic derivatives on continuous signals
before running a synchrony estimator.

- [`smooth_signal()`](https://jmgirard.github.io/bsync/reference/smooth_signal.md)
  : Smooth a Time Series Signal
- [`trim_edges()`](https://jmgirard.github.io/bsync/reference/trim_edges.md)
  : Trim Edge Effects from Data
- [`calc_velocity_1d()`](https://jmgirard.github.io/bsync/reference/calc_velocity_1d.md)
  : Calculate 1D Velocity
- [`calc_speed_1d()`](https://jmgirard.github.io/bsync/reference/calc_speed_1d.md)
  : Calculate 1D Speed
- [`calc_speed_2d()`](https://jmgirard.github.io/bsync/reference/calc_speed_2d.md)
  : Calculate 2D Speed
- [`calc_speed_3d()`](https://jmgirard.github.io/bsync/reference/calc_speed_3d.md)
  : Calculate 3D Speed
- [`evaluate_signal_power()`](https://jmgirard.github.io/bsync/reference/evaluate_signal_power.md)
  : Evaluate Signal Power and Suggest Downsampling Rate
- [`plot(`*`<signal_power_res>`*`)`](https://jmgirard.github.io/bsync/reference/plot.signal_power_res.md)
  : Plot method for signal_power_res objects

## Preprocessing — resampling & gaps

Resample continuous signals and diagnose or repair time-series gaps.

- [`downsample_signal()`](https://jmgirard.github.io/bsync/reference/downsample_signal.md)
  : Downsample a Time Series Signal via Rolling Aggregation
- [`aggregate_by_time()`](https://jmgirard.github.io/bsync/reference/aggregate_by_time.md)
  : Aggregate Time Series Data by Time Bins
- [`diagnose_ts_gaps()`](https://jmgirard.github.io/bsync/reference/diagnose_ts_gaps.md)
  : Diagnose Missing Data Gaps in a Time Series
- [`impute_ts_gaps()`](https://jmgirard.github.io/bsync/reference/impute_ts_gaps.md)
  : Impute Missing Values in Continuous Time Series with Metadata

## Tidy interface

[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html), and
[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
methods for `bsync_surface` and `bsync_multiverse` objects, enabling
integration with tidyverse workflows.

- [`reexports`](https://jmgirard.github.io/bsync/reference/reexports.md)
  [`tidy`](https://jmgirard.github.io/bsync/reference/reexports.md)
  [`glance`](https://jmgirard.github.io/bsync/reference/reexports.md)
  [`as_tibble`](https://jmgirard.github.io/bsync/reference/reexports.md)
  : Objects exported from other packages
- [`tidy(`*`<bsync_surface>`*`)`](https://jmgirard.github.io/bsync/reference/tidy.bsync_surface.md)
  : Tidy a bsync_surface object into a tibble of per-cell results
- [`tidy(`*`<bsync_multiverse>`*`)`](https://jmgirard.github.io/bsync/reference/tidy.bsync_multiverse.md)
  : Tidy a bsync_multiverse object into the specification grid
- [`glance(`*`<bsync_surface>`*`)`](https://jmgirard.github.io/bsync/reference/glance.bsync_surface.md)
  : One-row summary of a bsync_surface object
- [`glance(`*`<bsync_multiverse>`*`)`](https://jmgirard.github.io/bsync/reference/glance.bsync_multiverse.md)
  : One-row robustness summary of a bsync_multiverse object
- [`as_tibble(`*`<bsync_surface>`*`)`](https://jmgirard.github.io/bsync/reference/as_tibble.bsync_surface.md)
  : Convert a bsync_surface object to a tibble
- [`as_tibble(`*`<bsync_multiverse>`*`)`](https://jmgirard.github.io/bsync/reference/as_tibble.bsync_multiverse.md)
  : Convert a bsync_multiverse object to a tibble

## Data

- [`sim_dyad`](https://jmgirard.github.io/bsync/reference/sim_dyad.md) :
  Simulated Dyadic 3D Positional Data

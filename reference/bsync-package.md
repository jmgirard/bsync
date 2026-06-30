# bsync: Behavioral Synchrony Analyses

A comprehensive toolkit for analyzing interpersonal and behavioral
synchrony in continuous dyadic time series. Provides efficient, C++
backed implementations of three windowed estimators: cross-correlation
(WCC), dynamic time warping (WDTW), and Granger causality, each with a
full surrogate (pseudo-synchrony) significance testing layer using
circular-shift and phase-randomization generators. Includes peak and
valley optima picking, a rolling leadership asymmetry index, and a
preprocessing pipeline featuring zero-phase smoothing, time-bin
aggregation, gap imputation, and kinematic velocity and speed
calculations. Theory-driven and data-driven hyperparameter helpers guide
window size and lag selection.

## See also

Useful links:

- <https://jmgirard.github.io/bsync/>

- <https://github.com/jmgirard/bsync>

- Report bugs at <https://github.com/jmgirard/bsync/issues>

## Author

**Maintainer**: Jeffrey M. Girard <me@jmgirard.com>
([ORCID](https://orcid.org/0000-0002-7359-3746))

Authors:

- Jeffrey M. Girard <me@jmgirard.com>
  ([ORCID](https://orcid.org/0000-0002-7359-3746))

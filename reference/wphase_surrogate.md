# Calculate Surrogate Windowed Phase Synchrony

Calculate Surrogate Windowed Phase Synchrony

## Usage

``` r
wphase_surrogate(
  x,
  y,
  y_surrogates,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1
)
```

## Arguments

- x:

  A numeric vector containing a time series.

- y:

  A numeric vector containing a time series.

- y_surrogates:

  A matrix of surrogate time series for `y` (columns are surrogates).

- time:

  An optional numeric vector representing the timestamps for the data.
  Default is `NULL`.

- window_size:

  A positive integer indicating the size of each window.

- lag_max:

  A positive integer indicating the maximum lag to try.

- window_increment:

  A positive integer indicating the window shift increment. Default is
  1.

- lag_increment:

  A positive integer indicating the lag shift increment. Default is 1.

## Value

A list object of class "wphase_surr".

## Details

The p-value is the add-one form p = (b + 1) / (n + 1) (Phipson & Smyth,
2010), where n is the number of surrogate columns and b is the number of
surrogates whose aggregate statistic is **at least as large as** the
observed statistic. The `+ 1` counts the observed data as one more draw
from the null, so the p-value is never 0. Its smallest value is 1 / (n +
1). The aggregate – the mean phase-locking value over all window x lag
combinations (`mean_plv`) – is computed identically on the observed data
and every surrogate via the same internal helper, so the null
distribution and the observed value are directly comparable (Invariant
2: surrogate nulls match the observed statistic).

Phase extraction (the analytic signal via
[`gsignal::hilbert()`](https://rdrr.io/pkg/gsignal/man/hilbert.html)) is
applied to each surrogate column exactly as to the observed series.
Circular-shift surrogates
([`generate_surrogate_circular()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_circular.md))
preserve each series' own phase continuity while breaking cross-series
alignment, which makes them a natural null for a phase statistic;
phase-randomization surrogates
([`generate_surrogate_phase()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_phase.md))
instead destroy the surrogate's phase structure itself – state which
null you are testing.

## References

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero: calculating exact p-values when permutations are randomly drawn.
*Statistical Applications in Genetics and Molecular Biology*, 9(1),
Article 39.
[doi:10.2202/1544-6115.1585](https://doi.org/10.2202/1544-6115.1585)

## Examples

``` r
# \donttest{
# Two-step pipeline: generate a null matrix, then test the observed synchrony
y_surr <- generate_surrogate_circular(sim_dyad$z_B, n_surrogates = 100)
res <- wphase_surrogate(
  x = sim_dyad$z_A,
  y = sim_dyad$z_B,
  y_surrogates = y_surr,
  window_size = 96,
  lag_max = 10
)
res
#> 
#> ── Windowed Phase Synchrony Surrogate Analysis (Pseudo-Synchrony) ──────────────
#> Permutations: 100
#> Observed Mean PLV: 0.9993
#> Average Null Mean PLV: 0.9791
#> Empirical p-value: 0.0099
#> ✔ Observed phase synchrony is significantly greater than chance.
#> ℹ Note: 100 permutations may be too few for stable p-values.
#> Consider setting `n_surrogates >= 1000` for final reporting.
# }
```

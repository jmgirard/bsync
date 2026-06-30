# Calculate Surrogate Windowed Dynamic Time Warping (WDTW)

Calculate Surrogate Windowed Dynamic Time Warping (WDTW)

## Usage

``` r
wdtw_surrogate(
  x,
  y,
  y_surrogates,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1,
  scale_method = c("global", "local", "none"),
  distance_metric = c("L2", "L1"),
  fast_method = FALSE
)
```

## Arguments

- x:

  A numeric vector containing the reference time series.

- y:

  A numeric vector containing the query time series.

- y_surrogates:

  A matrix of surrogate time series for \`y\` (columns are surrogates).

- time:

  An optional numeric vector representing the timestamps for the data.
  Default is \`NULL\`.

- window_size:

  A positive integer indicating the size of the rolling window.

- lag_max:

  A positive integer indicating the maximum lag to try.

- window_increment:

  A positive integer indicating the step size for the rolling window.
  Default is 1.

- lag_increment:

  A positive integer indicating the lag shift increment. Default is 1.

- scale_method:

  Character string specifying how to standardize the data. Default is
  \`"global"\`.

- distance_metric:

  Character string specifying the local cost function. Default is
  \`"L2"\`.

- fast_method:

  Logical. If \`TRUE\`, severely reduces computation time by only
  evaluating surrogate alignments at lag 0. \*\*See Details for the
  statistical caveat.\*\* Default is \`FALSE\`.

## Value

A list object of class "wdtw_surr".

## Details

The p-value is the proportion of surrogates whose aggregate statistic is
\*\*at most as large as\*\* the observed statistic (lower DTW distance =
better alignment). The aggregate is \`mean(dtw_dist)\` over all window ×
lag combinations — the same quantity stored in
\`wdtw_res\$aggregate\[\["mean_distance"\]\]\` — computed identically on
both the observed data and every surrogate, so the null distribution and
the observed value are directly comparable.

\*\*\`fast_method\` warning:\*\* when \`fast_method = TRUE\`, surrogates
are evaluated at lag 0 only, while the observed statistic is computed
over all lags. The null and observed aggregates therefore cover
different lag ranges, making the resulting p-value approximate. Use only
for quick exploratory checks, never for reporting.

## Examples

``` r
# \donttest{
# DTW runs n_surrogates + 1 times, so this example uses a short subset and a
# small surrogate count for speed; use the full series and
# n_surrogates >= 1000 for reporting.
xs <- sim_dyad$x_A[1:250]
ys <- sim_dyad$x_B[1:250]
y_surr <- generate_surrogate_circular(ys, n_surrogates = 19)
res <- wdtw_surrogate(
  x = xs,
  y = ys,
  y_surrogates = y_surr,
  window_size = 50,
  lag_max = 5
)
res
#> 
#> ── WDTW Surrogate Analysis (Pseudo-Synchrony) ──────────────────────────────────
#> Permutations: 19
#> Observed Mean Cost: 27.8222
#> Average Null Cost: 30.1765
#> Empirical p-value: < 0.0526315789473684
#> ✔ Observed cost is significantly lower than chance (stronger alignment).
#> ℹ Note: 19 permutations may be too few for stable p-values.
#> Consider setting `n_surrogates >= 1000` for final reporting.
# }
```

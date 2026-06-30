# Calculate Surrogate Windowed Cross-Correlations

Calculate Surrogate Windowed Cross-Correlations

## Usage

``` r
wcc_surrogate(
  x,
  y,
  y_surrogates,
  time = NULL,
  window_size,
  lag_max,
  window_increment = 1,
  lag_increment = 1,
  na.rm = TRUE,
  statistic = c("mean_abs_z", "peak")
)
```

## Arguments

- x:

  A numeric vector containing a time series.

- y:

  A numeric vector containing a time series.

- y_surrogates:

  A matrix of surrogate time series for \`y\` (columns are surrogates).

- time:

  An optional numeric vector representing the timestamps for the data.
  Default is \`NULL\`.

- window_size:

  A positive integer indicating the size of each window.

- lag_max:

  A positive integer indicating the maximum lag to try.

- window_increment:

  A positive integer indicating the window shift increment. Default is
  1.

- lag_increment:

  A positive integer indicating the lag shift increment. Default is 1.

- na.rm:

  A logical indicating whether to remove missing values. Default is
  \`TRUE\`.

- statistic:

  A character string specifying the aggregate statistic; must match the
  value passed to \`wcc()\`. \`"mean_abs_z"\` (default) or \`"peak"\`.
  See \`wcc()\` for details.

## Value

A list object of class "wcc_surr".

## Details

The p-value is the proportion of surrogates whose aggregate statistic is
\*\*at least as large as\*\* the observed statistic. The aggregate —
either \`"mean_abs_z"\` or \`"peak"\` — is computed identically on the
observed data and every surrogate via the same internal helper, so the
null distribution and the observed value are guaranteed to be directly
comparable (Invariant 2: surrogate nulls match the observed statistic).

Pass the same \`statistic\` value you used in \`wcc()\` so that
\`observed_z\` and the surrogate draws use the same quantity.

## Examples

``` r
# \donttest{
# Two-step pipeline: generate a null matrix, then test the observed WCC
y_surr <- generate_surrogate_circular(sim_dyad$x_B, n_surrogates = 100)
res <- wcc_surrogate(
  x = sim_dyad$x_A,
  y = sim_dyad$x_B,
  y_surrogates = y_surr,
  window_size = 96,
  lag_max = 10
)
res
#> 
#> ── WCC Surrogate Analysis (Pseudo-Synchrony) ───────────────────────────────────
#> Permutations: 100
#> Observed Mean Abs. Fisher's Z: 0.09
#> Average Null Mean Abs. Fisher's Z: 0.0826
#> Empirical p-value: < 0.01
#> ✔ Observed synchrony is significantly greater than chance.
#> ℹ Note: 100 permutations may be too few for stable p-values.
#> Consider setting `n_surrogates >= 1000` for final reporting.
# }
```

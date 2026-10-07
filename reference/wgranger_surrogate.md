# Calculate Surrogate Windowed Granger Causality

Calculate Surrogate Windowed Granger Causality

## Usage

``` r
wgranger_surrogate(
  x,
  y,
  y_surrogates,
  time = NULL,
  window_size,
  ar_order = 1,
  window_increment = 1
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

  A positive integer indicating the size of the rolling window.

- ar_order:

  A positive integer specifying the Autoregressive (AR) order. Default
  is 1.

- window_increment:

  A positive integer indicating the step size for the rolling window.
  Default is 1.

## Value

A list object of class "wgranger_surr".

## Details

Two p-values are returned: one for x → y and one for y → x. Each is the
add-one form p = (b + 1) / (n + 1) (Phipson & Smyth, 2010), where n is
the number of surrogate columns and b is the number of surrogates whose
mean F-statistic across windows is **at least as large as** the
corresponding observed mean F-statistic. The `+ 1` counts the observed
data as one more draw from the null, so neither p-value is ever 0. The
smallest value is 1 / (n + 1). Both null distributions are built with
the same aggregate (`mean(f_xy)` and `mean(f_yx)`) as the observed
statistics stored in `wgranger_res$results_df`, so the null and observed
values are directly comparable.

## References

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero: calculating exact p-values when permutations are randomly drawn.
*Statistical Applications in Genetics and Molecular Biology*, 9(1),
Article 39.
[doi:10.2202/1544-6115.1585](https://doi.org/10.2202/1544-6115.1585)

## Examples

``` r
# \donttest{
y_surr <- generate_surrogate_circular(sim_dyad$x_B, n_surrogates = 100)
res <- wgranger_surrogate(
  x = sim_dyad$x_A,
  y = sim_dyad$x_B,
  y_surrogates = y_surr,
  window_size = 96
)
res
#> 
#> ── Windowed Granger Surrogate Analysis ─────────────────────────────────────────
#> 
#> ── Direction: x -> y ──
#> 
#> Permutations: 100
#> Observed Mean F-statistic: 1.0923
#> Average Null F-statistic: 1.0122
#> Empirical p-value: 0.3366
#> ! Predictive power (x -> y) is not significantly different from chance.
#> 
#> ── Direction: y -> x ──
#> 
#> Observed Mean F-statistic: 1.8485
#> Average Null F-statistic: 1.0024
#> Empirical p-value: 0.0099
#> ✔ Predictive power (y -> x) is significantly greater than chance.
#> ℹ Note: 100 permutations may be too few for stable p-values.
#> Consider setting `n_surrogates >= 1000` for final reporting.
# }
```

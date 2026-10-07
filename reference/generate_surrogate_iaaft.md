# Generate IAAFT Surrogates

Builds surrogates with the iterative amplitude-adjusted Fourier
transform (IAAFT) of Schreiber and Schmitz (1996). Each surrogate is a
new series that matches both the power spectrum and the value
distribution of `y`. One of the two matches exactly and the other
closely, as `match` selects.

## Usage

``` r
generate_surrogate_iaaft(
  y,
  n_surrogates = 100,
  max_iter = 1000,
  match = c("spectrum", "values")
)
```

## Arguments

- y:

  A numeric vector with at least 3 finite values and no `NA`.

- n_surrogates:

  A single positive integer: the number of surrogates. Default is `100`.

- max_iter:

  A single positive integer: the maximum number of iterations per
  surrogate. Default is `1000`.

- match:

  Which property matches exactly: `"spectrum"` (default), the Fourier
  amplitudes of `y`, or `"values"`, the values of `y`. The other
  property matches closely. See Details.

## Value

A numeric matrix with `length(y)` rows and `n_surrogates` columns, ready
to pass as `y_surrogates` to
[`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md)
and the other surrogate wrappers.

## Details

**Which null this tests.** Phase randomization
([`generate_surrogate_phase()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_phase.md))
keeps the power spectrum but makes the values Gaussian, so a skewed or
bounded signal gets surrogates with a different value distribution.
IAAFT keeps the spectrum and the value distribution. The null is a
Gaussian linear process seen through a fixed, monotone transform, for
example a skewed movement-energy signal. Against this null, a
significant result cannot come from the shape of the value distribution
or from the autocorrelation of `y` alone.

**Algorithm.** Each surrogate starts from a random shuffle of `y`. Each
iteration takes the Fourier transform, replaces its amplitudes with
those of `y` while it keeps the phases, and transforms back. This is the
spectrum-adjusted series. The iteration then rank-orders it, so that the
series takes exactly the values of `y`. A surrogate stops when the
rank-ordering no longer changes it, which is where the iteration ends
(Schreiber & Schmitz, 1996, p. 2). Apart from a cyclic shift of `y`, a
finite series in general cannot match the spectrum and the values both
exactly (p. 2). A very short series has few orderings, and its
surrogates can be cyclic shifts or reversals of `y`, which break little
of its structure.

**`match`.** `"spectrum"` (the default) returns the last
spectrum-adjusted series. Its Fourier amplitudes equal those of `y`
exactly, and its values are close to the values of `y`: closer than the
values of a phase surrogate, in the package tests on a skewed
autoregressive series. `"values"` returns the rank-ordered series, as in
the paper. It holds exactly the values of `y`, and its spectrum is close
to that of `y`. The rank step lowers the autocorrelation slightly, and a
synchrony statistic that grows with autocorrelation, such as the mean
absolute Fisher z of
[`wcc()`](https://jmgirard.github.io/bsync/reference/wcc.md), then reads
high against these surrogates. So a test against `"values"` surrogates
can reject a true null more often than its nominal level, and the
spectrum form is the default. In a size check run on 2026-10-06
(independent pairs of AR(1) series x_n = 0.7 x\_(n-1) + e_n observed as
x_n^3, length 512, 19 surrogates per pair,
[`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md)
with `window_size = 32`, `lag_max = 4`, `window_increment = 16`), the
spectrum form rejected 17 of 400 pairs (4.25%) at p \<= .05. On the
first 200 of those pairs, the values form rejected 21 (10.5%) and the
spectrum form 9 (4.5%). The package's own tests repeat this check.

**`max_iter`.** The spectral error falls about as 1 / i over the first
iterations (p. 3), and the paper's Fig. 2 runs to 1000 iterations (p.
2). The default `1000` is a cap on run time, not a target. A surrogate
that reaches `max_iter` without a fixed point is still returned, and the
call warns with the number of such surrogates. Long series need more
iterations: on one AR(1)-cube series of 10000 samples (run 2026-10-06),
4 of 5 surrogates had not reached the fixed point at 1000 iterations.
With the default `match`, such surrogates still match the spectrum
exactly.

**Missing values.** The Fourier transform cannot take `NA`, so `y` must
be complete. Fill gaps first, for example with
[`impute_ts_gaps()`](https://jmgirard.github.io/bsync/reference/impute_ts_gaps.md).

## References

Schreiber, T., & Schmitz, A. (1996). Improved surrogate data for
nonlinearity tests. *Physical Review Letters*, 77(4), 635-638.
[doi:10.1103/PhysRevLett.77.635](https://doi.org/10.1103/PhysRevLett.77.635)

## See also

[`generate_surrogate_phase()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_phase.md)
and
[`generate_surrogate_circular()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_circular.md)
for the other within-dyad nulls;
[`generate_surrogate_pseudo()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_pseudo.md)
for the between-dyad null;
[`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md),
[`wdtw_surrogate()`](https://jmgirard.github.io/bsync/reference/wdtw_surrogate.md),
[`wgranger_surrogate()`](https://jmgirard.github.io/bsync/reference/wgranger_surrogate.md),
[`wphase_surrogate()`](https://jmgirard.github.io/bsync/reference/wphase_surrogate.md).

## Examples

``` r
# Build 100 IAAFT surrogates: exact spectrum, close values
surr <- generate_surrogate_iaaft(sim_dyad$x_B, n_surrogates = 100)
dim(surr)
#> [1] 2400  100
all.equal(Mod(fft(surr[, 1])), Mod(fft(sim_dyad$x_B)))
#> [1] TRUE

# The "values" form: exact values, close spectrum
surr_v <- generate_surrogate_iaaft(sim_dyad$x_B, 10, match = "values")
all.equal(sort(surr_v[, 1]), sort(sim_dyad$x_B))
#> [1] TRUE
```

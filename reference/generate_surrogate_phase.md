# Generate Phase-Randomized Surrogates (Fourier Transform)

Generate Phase-Randomized Surrogates (Fourier Transform)

## Usage

``` r
generate_surrogate_phase(y, n_surrogates = 100, trim_odd = FALSE)
```

## Arguments

- y:

  A numeric vector containing a time series.

- n_surrogates:

  Integer specifying the number of surrogates. Default is 100.

- trim_odd:

  Logical. Odd-length series are supported natively; when \`TRUE\`, the
  final observation is dropped instead (legacy behavior). Default is
  \`FALSE\`.

## Value

A matrix where each column is a surrogate time series.

## Details

Randomizes the phases of the Fourier transform of \`y\` while preserving
its amplitude spectrum (and therefore its autocovariance / power
spectrum). The DC term and – for even-length series – the Nyquist term
are real-valued components carrying a sign; they are preserved exactly
(not just their modulus), so the surrogate mean and variance match the
original for any signal, including negative-mean signals. Conjugate
symmetry is enforced on the remaining bins so the inverse transform is
real. Both even- and odd-length series are handled natively.

## See also

[`generate_surrogate_iaaft`](https://jmgirard.github.io/bsync/reference/generate_surrogate_iaaft.md),
which also keeps the value distribution, and
[`generate_surrogate_circular`](https://jmgirard.github.io/bsync/reference/generate_surrogate_circular.md)
for the other within-dyad nulls;
[`generate_surrogate_pseudo`](https://jmgirard.github.io/bsync/reference/generate_surrogate_pseudo.md)
for the between-dyad null.

## Examples

``` r
# Build 100 phase-randomized surrogates (preserves the power spectrum)
surr <- generate_surrogate_phase(sim_dyad$x_B, n_surrogates = 100)
dim(surr)
#> [1] 2400  100
```

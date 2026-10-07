# Generate Segment-Shuffling Surrogates

Cuts `y` into segments of `segment_size` samples and builds each
surrogate by putting the segments in a new order. Each segment keeps its
own values, autocorrelation, and any trend inside it. Only the time at
which each segment occurred changes.

## Usage

``` r
generate_surrogate_segment(y, segment_size, n_surrogates = 100)
```

## Arguments

- y:

  A numeric vector containing a time series. `NA` values are allowed.

- segment_size:

  A single whole number from 2 to `floor(length(y) / 2)`: the segment
  length in samples.

- n_surrogates:

  A single positive integer: the number of surrogates. Default is `100`.
  Above 8 segments, it must be below k! - 1.

## Value

A numeric matrix with `length(y)` rows and one column per surrogate,
ready to pass as `y_surrogates` to
[`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md)
and the other surrogate wrappers. It has `n_surrogates` columns, or k! -
1 columns if that is fewer.

## Details

**Which null this tests.** Each window of `x` has one segment of `y`
that occurred alongside it. The test asks whether that segment resembles
the window more than the other segments of `y` do. Under the null, the
time at which each segment of `y` occurred carries no information about
`x`, so any order of the segments is as good as the real one. The test
is exact if the segments of `y` are exchangeable, which is close to true
when a segment is long relative to the memory of `y`, and if every
window lies inside one segment (see "Windows and segments").

**Segments and orders.** The series has k =
`floor(length(y) / segment_size)` segments, cut from sample 1. Segment i
covers samples (i - 1) \* `segment_size` + 1 to i \* `segment_size`.
Each surrogate takes an order drawn uniformly from all k! - 1 orders
other than the original, and no order is drawn twice. So with 3 or more
segments, a surrogate can leave a segment in place: on average (k! - k)
/ (k! - 1) of the k segments, which is close to one from 4 segments on.
This rule makes the add-one p-value of the surrogate wrappers exact.
From 4 segments on, orders that move every segment do not form a group
with the original order, and a test against them rejects a true null too
often. With up to 8 segments and an `n_surrogates` of at least k! - 1,
the call returns all k! - 1 orders and says so. With 3 segments or
fewer, k! - 1 is below 19, so no test can reach p \<= .05, and the call
warns.

**The tail.** The last `length(y) - k * segment_size` samples do not
fill a segment. They stay in place at the end of every surrogate, and
the call says how many there are. Missing values move with their
segments.

**Windows and segments.** A window that crosses a segment boundary is
continuous in `y` but not in a surrogate, and the test then rejects a
true null too often. Every window and every lagged window lies inside
one segment under these settings:

- For [`wcc()`](https://jmgirard.github.io/bsync/reference/wcc.md),
  [`wdtw()`](https://jmgirard.github.io/bsync/reference/wdtw.md), and
  [`wphase()`](https://jmgirard.github.io/bsync/reference/wphase.md):
  `segment_size >= window_size + 2 * lag_max` and
  `window_increment = segment_size`.

- For
  [`wgranger()`](https://jmgirard.github.io/bsync/reference/wgranger.md):
  `segment_size = window_size = window_increment`.

The window grid then has k - 1 windows, or k windows when the tail is at
least `window_size + 2 * lag_max - 1` samples (`window_size - 1` for
[`wgranger()`](https://jmgirard.github.io/bsync/reference/wgranger.md)).
With `segment_size = window_size + 2 * lag_max`, that is a tail of
`segment_size - 1` samples. With k - 1 windows, the last segment serves
only as a source for the surrogates.
[`wphase()`](https://jmgirard.github.io/bsync/reference/wphase.md) takes
the Hilbert transform of the whole series, so a segment test of phase
synchrony is approximate even with aligned windows.
[`synchrony_multiverse()`](https://jmgirard.github.io/bsync/reference/synchrony_multiverse.md)
and
[`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md)
use these settings for `surrogate_method = "segment"`.

**SUSY.** The segment cut follows the `susy()` function of the SUSY
package (Tschacher & Meier, 2020). SUSY correlates each segment of one
series with each segment of the other and compares the real pairs with
the pairs of different segments. bsync instead builds full-length
reordered series, so that every surrogate wrapper can use them. SUSY's
pseudo pairs never match a segment with itself, but bsync's orders can
leave a segment in place, which keeps the p-value exact. SUSY drops the
tail, and bsync keeps it. bsync does not reproduce SUSY's effect size or
its numbers. SUSY's rule `segment >= 2 * maxlag` has the analogue
`segment_size >= window_size + 2 * lag_max`.

## References

Tschacher, W., & Meier, D. (2020). Physiological synchrony in
psychotherapy sessions. *Psychotherapy Research*, 30(5), 558-573.
[doi:10.1080/10503307.2019.1612114](https://doi.org/10.1080/10503307.2019.1612114)

## See also

[`generate_surrogate_circular()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_circular.md),
[`generate_surrogate_phase()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_phase.md),
and
[`generate_surrogate_iaaft()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_iaaft.md)
for the other within-dyad nulls;
[`generate_surrogate_pseudo()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_pseudo.md)
for the between-dyad null;
[`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md),
[`wdtw_surrogate()`](https://jmgirard.github.io/bsync/reference/wdtw_surrogate.md),
[`wgranger_surrogate()`](https://jmgirard.github.io/bsync/reference/wgranger_surrogate.md),
[`wphase_surrogate()`](https://jmgirard.github.io/bsync/reference/wphase_surrogate.md).

## Examples

``` r
# Windows of 96 samples with lags up to 10 need segments of at least
# 96 + 2 * 10 = 116 samples: 20 segments of sim_dyad's 2400 samples.
surr <- generate_surrogate_segment(
  sim_dyad$x_B,
  segment_size = 116,
  n_surrogates = 19
)
#> The last 80 samples of `y` do not fill a whole segment and stay in place.
dim(surr)
#> [1] 2400   19

# \donttest{
# Step the windows one segment at a time, so each lies inside a segment
wcc_surrogate(
  x = sim_dyad$x_A, y = sim_dyad$x_B, y_surrogates = surr,
  window_size = 96, lag_max = 10, window_increment = 116
)
#> 
#> ── WCC Surrogate Analysis (Pseudo-Synchrony) ───────────────────────────────────
#> Permutations: 19
#> Observed Mean Abs. Fisher's Z: 0.0868
#> Average Null Mean Abs. Fisher's Z: 0.0842
#> Empirical p-value: 0.4
#> ! Observed synchrony is not significantly different from chance.
#> ℹ Note: 19 permutations may be too few for stable p-values.
#> Consider setting `n_surrogates >= 1000` for final reporting.
# }
```

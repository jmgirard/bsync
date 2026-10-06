# Generate Pseudo-Dyad Surrogates for One Dyad

Builds a surrogate matrix for one target dyad whose columns are partner
series taken from *other* dyads in the sample. Pairing the target's `x`
with a partner who never interacted with that person gives a
pseudo-synchrony null: it keeps any co-movement that the shared task or
setting produces at the same moment, and removes the coupling that is
specific to the real interaction.

## Usage

``` r
generate_surrogate_pseudo(
  dyad_list,
  dyad,
  n_surrogates = NULL,
  keep_roles = TRUE
)
```

## Arguments

- dyad_list:

  A list with one element per dyad, the form
  [`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md)
  also reads. Each element is a data frame or a list. If its names
  contain `x` and `y` exactly once each, those two elements are read by
  name, in any position. If its names contain neither `x` nor `y`, or it
  has no names, its first two elements are read as `x` and `y`. Any
  other use of the names `x` and `y` is an error. Names are never
  matched partially. If `dyad_list` itself is named, the names identify
  the dyads in the output, so every name must be non-empty, not `NA`,
  and unique.

- dyad:

  A single integer: the index of the target dyad in `dyad_list`.

- n_surrogates:

  `NULL` (default) to use every eligible partner, or a single positive
  integer to draw that many without replacement.

- keep_roles:

  Logical. `TRUE` (default) draws partners only from the `y` series of
  other dyads. `FALSE` also allows their `x` series.

## Value

A numeric matrix with `length(y)` rows and one column per partner, ready
to pass as `y_surrogates` to
[`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md)
and the other surrogate wrappers. The attribute `"sources"` is a data
frame with one row per column and the columns `dyad` (the source dyad's
index in `dyad_list`), `role` (`"x"` or `"y"`), and `name` (the source
dyad's name in `dyad_list`, or `NA` if `dyad_list` has no names). If
`dyad_list` is named, the matrix column names are `"<name>:<role>"`.
Subsetting the matrix drops the attribute, so read it before you subset.

## Details

**Which null this tests.** Circular-shift and phase surrogates break all
time alignment, so a shared task structure (for example, both people
move when a trial starts) can look like synchrony against them. A
pseudo-dyad partner was recorded under the same task timeline but in a
different interaction, so the null keeps task-driven co-movement. A
significant result then means more synchrony than people in the same
setting show with a stranger (Kleinbub & Ramseyer, 2020).

**Length handling.** Partners are cropped *start-aligned*: each column
is the first `length(y)` samples of its source, so sample `t` of the
partner sits at the same task time as sample `t` of the target (the rMEA
convention). Partners shorter than the target's `y` cannot fill the
matrix and are excluded with a warning. How many longer partners were
cropped is announced. Missing values are passed through unchanged. Start
alignment assumes that every dyad was sampled at the same rate and that
sample 1 of every series is the same task onset. Resample and trim first
if not; this function cannot check it.

**Roles.** With `keep_roles = TRUE` (the default) every partner is the
`y` of another dyad. bsync's lead-lag sign depends on which series is
`x`, so keeping roles keeps the null comparable to the observed dyad.
This differs from rMEA's `shuffle()`, which mixes roles by default. Use
`keep_roles = FALSE` when partners are interchangeable: the `x` and the
`y` of every other dyad are then candidates, which doubles the pool.

**How many surrogates.** A sample of N dyads offers at most N - 1
partners per target (2(N - 1) with `keep_roles = FALSE`). The wrappers
report the add-one p-value (b + 1) / (n + 1), where n is the number of
partner columns and b counts the partners that score at least as high as
the real partner. The smallest possible p-value is 1 / (n + 1). With N
dyads and `keep_roles = TRUE`, n is at most N - 1, so the single-dyad
p-value is never below 1 / N, and `p <= .05` needs at least 20 dyads.
With `keep_roles = FALSE`, n is at most 2(N - 1), so the floor is 1 /
(2N - 1) and `p <= .05` needs at least 11 dyads. When every partner is
used with `keep_roles = TRUE`, the p-value is the real partner's rank
among the N series divided by N. If the real partner is no different
from a stranger, each rank has chance 1 / N, so the test keeps its
nominal size. With too few dyads to reach `p <= .05`, prefer the
sample-level comparison of
[`generate_pseudo_dyads()`](https://jmgirard.github.io/bsync/reference/generate_pseudo_dyads.md).
`n_surrogates = NULL` (the default) uses every eligible partner, with no
random draw. An integer draws that many partners without replacement.
Asking for more partners than exist is an error, because repeated
partners add no information to the null.

## References

Kleinbub, J. R., & Ramseyer, F. T. (2020). rMEA: An R package to assess
nonverbal synchronization in motion energy analysis time-series.
*Psychotherapy Research*.
[doi:10.1080/10503307.2020.1844334](https://doi.org/10.1080/10503307.2020.1844334)

## See also

[`generate_pseudo_dyads()`](https://jmgirard.github.io/bsync/reference/generate_pseudo_dyads.md)
for the sample-wide set of pseudo-dyads;
[`generate_surrogate_circular()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_circular.md)
and
[`generate_surrogate_phase()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_phase.md)
for within-dyad nulls;
[`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md),
[`wdtw_surrogate()`](https://jmgirard.github.io/bsync/reference/wdtw_surrogate.md),
[`wgranger_surrogate()`](https://jmgirard.github.io/bsync/reference/wgranger_surrogate.md),
[`wphase_surrogate()`](https://jmgirard.github.io/bsync/reference/wphase_surrogate.md).

## Examples

``` r
# Three "dyads" built from sim_dyad's axes (a stand-in for a real sample)
dyads <- list(
  list(x = sim_dyad$x_A, y = sim_dyad$x_B),
  list(x = sim_dyad$y_A, y = sim_dyad$y_B),
  list(x = sim_dyad$z_A, y = sim_dyad$z_B)
)
y_pseudo <- generate_surrogate_pseudo(dyads, dyad = 3)
dim(y_pseudo)
#> [1] 2400    2
attr(y_pseudo, "sources")
#>   dyad role name
#> 1    1    y <NA>
#> 2    2    y <NA>

# \donttest{
# Test dyad 3 against partners from the other dyads
wcc_surrogate(
  x = sim_dyad$z_A, y = sim_dyad$z_B, y_surrogates = y_pseudo,
  window_size = 96, lag_max = 10
)
#> 
#> ── WCC Surrogate Analysis (Pseudo-Synchrony) ───────────────────────────────────
#> Permutations: 2
#> Observed Mean Abs. Fisher's Z: 1.3057
#> Average Null Mean Abs. Fisher's Z: 0.0871
#> Empirical p-value: 0.3333
#> ! Observed synchrony is not significantly different from chance.
#> ℹ With 2 surrogates, the smallest possible p-value is 0.333, so no result can reach p <= .05.
#> ℹ Note: 2 permutations may be too few for stable p-values.
#> Consider setting `n_surrogates >= 1000` for final reporting.
# }
```

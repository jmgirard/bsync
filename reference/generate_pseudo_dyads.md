# Generate the Sample-Wide Set of Pseudo-Dyads

Pairs series from *different* dyads to build pseudo-dyads: pairs of
people who never interacted. Computing the same synchrony statistic on
the real dyads and on these pseudo-dyads gives the classic
pseudo-synchrony comparison: real dyads should be more synchronous than
pseudo-dyads recorded under the same task (Kleinbub & Ramseyer, 2020).

## Usage

``` r
generate_pseudo_dyads(dyad_list, n_pairs = NULL, keep_roles = TRUE)
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

- n_pairs:

  `NULL` (default) to return every pairing, or a single positive integer
  to draw that many without replacement.

- keep_roles:

  Logical. `TRUE` (default) pairs the `x` of one dyad with the `y` of
  another. `FALSE` pairs any two series from different dyads.

## Value

A list of pseudo-dyads. Each element is `list(x = , y = )`, the dyad
form that
[`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md)
and
[`generate_surrogate_pseudo()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_pseudo.md)
read. The attribute `"sources"` is a data frame with one row per
pseudo-dyad and the columns `x_dyad`, `x_role`, `x_name`, `y_dyad`,
`y_role`, and `y_name`: the index, role, and name in `dyad_list` of each
series. The name columns are `NA` if `dyad_list` has no names. If
`dyad_list` is named, the list elements are named
`"<x_name>:<x_role>|<y_name>:<y_role>"`. These element names are labels:
if dyad names contain `:` or `|`, two labels can be equal, and the list
cannot then be passed back as a named `dyad_list`. Subsetting the list
drops the attribute, so read it before you subset.

## Details

**Pairings.** With `keep_roles = TRUE` (the default), each pseudo-dyad
is the `x` of dyad i with the `y` of dyad j, for every i != j: N(N - 1)
ordered pairs for N dyads. bsync's lead-lag sign depends on which series
is `x`, so keeping roles keeps pseudo-dyads comparable to real ones.
With `keep_roles = FALSE`, any two series from different dyads form a
pair, in either role: 2N(N - 1) unordered pairs. This is the pair set of
rMEA's `shuffle()`, which mixes roles by default, with the same
assignment of series to `x` and `y`: a pair of an `x` series and a `y`
series keeps the `x` series as `x`.

**Length handling.** Both series of a pseudo-dyad are cropped to the
shorter one, start-aligned (the rMEA convention), so sample `t` of each
series sits at the same task time. The number of cropped pseudo-dyads is
announced. Missing values are passed through unchanged. Start alignment
assumes that every dyad was sampled at the same rate and that sample 1
of every series is the same task onset. Resample and trim first if not;
this function cannot check it.

**How many.** `n_pairs = NULL` (the default) returns every pairing in a
fixed order, with no random draw. An integer draws that many pairings
without replacement. Asking for more than exist is an error.

## References

Kleinbub, J. R., & Ramseyer, F. T. (2020). rMEA: An R package to assess
nonverbal synchronization in motion energy analysis time-series.
*Psychotherapy Research*.
[doi:10.1080/10503307.2020.1844334](https://doi.org/10.1080/10503307.2020.1844334)

## See also

[`generate_surrogate_pseudo()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_pseudo.md)
for a per-dyad surrogate matrix, which
[`wcc_surrogate()`](https://jmgirard.github.io/bsync/reference/wcc_surrogate.md),
[`wdtw_surrogate()`](https://jmgirard.github.io/bsync/reference/wdtw_surrogate.md),
[`wgranger_surrogate()`](https://jmgirard.github.io/bsync/reference/wgranger_surrogate.md),
and
[`wphase_surrogate()`](https://jmgirard.github.io/bsync/reference/wphase_surrogate.md)
accept;
[`generate_surrogate_circular()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_circular.md),
[`generate_surrogate_phase()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_phase.md),
[`generate_surrogate_iaaft()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_iaaft.md),
and
[`generate_surrogate_segment()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_segment.md)
for within-dyad nulls;
[`wcc()`](https://jmgirard.github.io/bsync/reference/wcc.md),
[`wdtw()`](https://jmgirard.github.io/bsync/reference/wdtw.md),
[`wgranger()`](https://jmgirard.github.io/bsync/reference/wgranger.md),
and [`wphase()`](https://jmgirard.github.io/bsync/reference/wphase.md)
to compute the statistic on each pseudo-dyad.

## Examples

``` r
# Three "dyads" built from sim_dyad's axes (a stand-in for a real sample)
dyads <- list(
  list(x = sim_dyad$x_A, y = sim_dyad$x_B),
  list(x = sim_dyad$y_A, y = sim_dyad$y_B),
  list(x = sim_dyad$z_A, y = sim_dyad$z_B)
)
pseudo <- generate_pseudo_dyads(dyads)
length(pseudo) # 3 * 2 = 6 role-keeping pairings
#> [1] 6
attr(pseudo, "sources")
#>   x_dyad x_role x_name y_dyad y_role y_name
#> 1      1      x   <NA>      2      y   <NA>
#> 2      1      x   <NA>      3      y   <NA>
#> 3      2      x   <NA>      1      y   <NA>
#> 4      2      x   <NA>      3      y   <NA>
#> 5      3      x   <NA>      1      y   <NA>
#> 6      3      x   <NA>      2      y   <NA>

# Mean |Fisher z| on each pseudo-dyad: the pseudo-synchrony baseline
pseudo_z <- vapply(pseudo, function(d) {
  wcc(d$x, d$y, window_size = 96, lag_max = 10, window_increment = 48)$
    aggregate[[1]]
}, numeric(1))
summary(pseudo_z)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.07296 0.07693 0.08187 0.08081 0.08520 0.08657 
```

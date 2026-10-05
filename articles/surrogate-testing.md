# Surrogate Testing for Significance

In interpersonal synchrony analysis, finding a high correlation or a low
alignment cost is only the first step. Because behavioral time series
are inherently autocorrelated, standard statistical thresholds (like \\r
\> 0\\) are often misleading. A high correlation might occur purely by
chance because both signals share the same slow-moving trends.

Surrogate testing allows us to build an empirical null distribution—a
“chance baseline”—to determine if our observed results are truly
meaningful.

## 1. Choosing a Surrogate Method

The **bsync** package provides three methods for generating null data.
The first two work within one dyad. The third needs a sample of dyads.
Choosing the right one depends on the nature of your signal and on the
null hypothesis you want to test.

### 1.1 Circular Shift (`generate_surrogate_circular`)

The circular shift method takes a time series and shifts it by a random
interval, wrapping the end of the series back to the beginning.

- **How it works:** It destroys the temporal alignment between Person A
  and Person B while perfectly preserving the exact amplitude
  distribution and the power spectrum of the original signal.
- **When to use it:** This is our default recommendation for most dyadic
  interaction data. It is computationally efficient and effectively
  breaks the cross-signal dependency while keeping the “rhythm” of the
  individual intact.
- **Important Note:** If your data has strong global trends (e.g., a
  person slowly moving across the room over 60 seconds), a simple
  circular shift might not break that trend. Ensure you use the
  `lag_max` argument to ensure shifts are large enough to decouple the
  signals.

### 1.2 Phase Randomization (`generate_surrogate_phase`)

Phase randomization uses the Fourier Transform to decompose a signal
into its frequency components, randomizes the phases, and then
reconstructs the signal.

- **How it works:** It creates a signal that has the exact same power
  spectrum (frequency content) as the original, but a totally randomized
  temporal structure.
- **When to use it:** This is the “gold standard” for highly oscillatory
  data (e.g., rhythmic postural sway, physiological signals like heart
  rate or respiration).
- **Important Note:** Because this method creates an entirely new signal
  from frequency components, it is statistically more conservative than
  circular shifting. It is significantly more computationally expensive
  to generate than circular shifts.

### 1.3 Pseudo-Dyads (`generate_surrogate_pseudo`, `generate_pseudo_dyads`)

A pseudo-dyad pairs one person’s series with a partner from a
*different* dyad: two people who never interacted. This is the
pseudo-synchrony approach of the rMEA package (Kleinbub & Ramseyer,
2020).

- **How it works:** The partner series is real data, recorded under the
  same task timeline. Series of unequal length are cropped
  start-aligned, so sample `t` of both series sits at the same task
  time.
  [`generate_pseudo_dyads()`](https://jmgirard.github.io/bsync/reference/generate_pseudo_dyads.md)
  crops each pair to its shorter series.
  [`generate_surrogate_pseudo()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_pseudo.md)
  crops longer partners to the length of the target’s `y`, and excludes
  shorter partners with a warning. Start alignment assumes that every
  dyad was sampled at the same rate and that sample 1 of every series is
  the same task onset. Resample and trim first if that is not true.
- **Which null it tests:** Circular shifts and phase randomization break
  all time alignment. If the task itself makes both people move at the
  same moments (a trial starts, a video plays), that shared movement
  looks like synchrony against those nulls. A pseudo-dyad partner went
  through the same task, so the pseudo-dyad null keeps task-driven
  co-movement. A significant result then means “more synchronous than
  with a stranger in the same setting.”
- **When to use it:** When you have a sample of dyads that did the same
  task, and you want to separate interaction-specific synchrony from
  synchrony that the task produces.
- **Important Note:** A sample of N dyads gives each target at most N -
  1 partners (2(N - 1) with `keep_roles = FALSE`), so the p-value of one
  dyad has a resolution of 1 / (N - 1). The sample-level comparison in
  Section 3 uses all N(N - 1) pseudo-dyads.

## 2. Theoretical Pipeline

Regardless of the generation method, the logic of the **bsync**
surrogate pipeline remains consistent:

1.  **Decouple:** Generate a matrix of \\N\\ null time series that
    remove the synchrony while keeping the individual characteristics
    (or, for pseudo-dyads, that replace the partner with a person from
    another dyad).
2.  **Evaluate:** Run your chosen synchrony metric (WCC, WDTW, or WGC)
    against every column in that matrix.
3.  **Compare:** Calculate the empirical p-value as the proportion of
    surrogate results that are as extreme as, or more extreme than, your
    observed result.

``` r

library(bsync)

# 1. Generate 1000 surrogates using your preferred method
null_matrix <- generate_surrogate_circular(
  y = dyad_data$person_B,
  n_surrogates = 1000
)

# 2. Compare observed results to the null distribution
results <- wcc_surrogate(
  x = dyad_data$person_A,
  y = dyad_data$person_B,
  y_surrogates = null_matrix,
  window_size = 90,
  lag_max = 45
)

# 3. Interpret the empirical p-value
print(results)
```

## 3. Pseudo-Dyads: A Worked Example

The pseudo-dyad generators read a *list of dyads*. Each element is a
data frame (first two columns are `x` and `y`) or a list with elements
`x` and `y`, the same form that
[`autotune_wcc()`](https://jmgirard.github.io/bsync/reference/autotune_wcc.md)
reads. Here we simulate a sample of 10 dyads that did the same task.
Every person follows one shared task timeline, plus independent noise.
The dyads have no interaction-specific coupling at all.

``` r

library(bsync)
set.seed(2026)

n_dyads <- 10
n <- 1200
smooth <- function(v, k) as.numeric(stats::filter(v, rep(1 / k, k), circular = TRUE))

# One task timeline that every participant follows
task <- as.numeric(scale(smooth(cumsum(rnorm(n)), 40)))

task_only <- lapply(seq_len(n_dyads), function(i) {
  list(x = task + rnorm(n), y = task + rnorm(n))
})
```

We test each dyad twice: against circular-shift surrogates and against
pseudo-dyad surrogates.
[`generate_surrogate_pseudo()`](https://jmgirard.github.io/bsync/reference/generate_surrogate_pseudo.md)
returns a matrix whose columns are the `y` series of the other nine
dyads, ready for any surrogate wrapper.

``` r

test_dyad <- function(i, y_surrogates) {
  d <- task_only[[i]]
  wcc_surrogate(
    x = d$x, y = d$y, y_surrogates = y_surrogates,
    window_size = 100, lag_max = 20, window_increment = 25
  )$p_value
}

p_circular <- vapply(seq_len(n_dyads), function(i) {
  test_dyad(i, generate_surrogate_circular(task_only[[i]]$y, 99, lag_max = 20))
}, numeric(1))

p_pseudo <- vapply(seq_len(n_dyads), function(i) {
  test_dyad(i, generate_surrogate_pseudo(task_only, dyad = i))
}, numeric(1))

# Dyads called significant at .05, out of 10
c(circular = sum(p_circular < 0.05), pseudo = sum(p_pseudo < 0.05))
#> circular   pseudo 
#>        9        2
```

The circular-shift null breaks the shared task timeline, so it reads
task-driven co-movement as synchrony. The pseudo-dyad null keeps that
timeline, so it flags far fewer dyads. With nine partners per dyad,
`p_value` can only be 0, 1/9, 2/9, and so on, and `p < .05` means that
the real partner beat all nine strangers. If the real partner is no
different from a stranger, that happens with chance 1/10, so a
single-dyad pseudo-dyad test at `p < .05` has a false-positive rate of
about 1/N for N dyads. With small samples, prefer the sample-level
comparison below.

For a sample-level comparison,
[`generate_pseudo_dyads()`](https://jmgirard.github.io/bsync/reference/generate_pseudo_dyads.md)
builds every role-keeping pairing (the `x` of one dyad with the `y` of
another: 10 × 9 = 90 pseudo-dyads). Compute the same statistic on the
real dyads and on the pseudo-dyads, then compare the two distributions.
Here we add a real coupling to each dyad: person B follows person A by
10 samples.

``` r

coupled <- lapply(seq_len(n_dyads), function(i) {
  own <- smooth(rnorm(n), 5)
  list(
    x = task + own + rnorm(n, sd = 0.5),
    y = task + c(rep(0, 10), own[1:(n - 10)]) + rnorm(n, sd = 0.5)
  )
})
pseudo <- generate_pseudo_dyads(coupled)
length(pseudo)
#> [1] 90
head(attr(pseudo, "sources"))
#>   x_dyad x_role y_dyad y_role
#> 1      1      x      2      y
#> 2      1      x      3      y
#> 3      1      x      4      y
#> 4      1      x      5      y
#> 5      1      x      6      y
#> 6      1      x      7      y

mean_abs_z <- function(d) {
  wcc(d$x, d$y, window_size = 100, lag_max = 20, window_increment = 25)$aggregate[[1]]
}
real_z <- vapply(coupled, mean_abs_z, numeric(1))
pseudo_z <- vapply(pseudo, mean_abs_z, numeric(1))

summary(real_z)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.1835  0.1974  0.2011  0.2052  0.2097  0.2458
summary(pseudo_z)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.1430  0.1620  0.1681  0.1700  0.1774  0.1982
```

The real dyads score higher than the pseudo-dyads, because only real
partners share the 10-sample coupling. Both groups share the task
timeline. How you compare the two distributions (a t-test, a permutation
test, or a mixed model) is your choice. The generators only build the
pseudo-dyads.

## 4. Best Practices

- **Permutation Count:** For a stable p-value, we recommend at least
  1,000 permutations. Using fewer (e.g., 100) is fine for rapid
  exploratory analysis, but it will lack the granularity needed to
  report significance levels like \\p \< .001\\.
- **Power and Speed:** Surrogate analysis is “embarrassingly parallel.”
  For heavy analyses like WDTW, always use the `future` package to
  distribute these 1,000+ calculations across your machine’s CPU cores.
- **Interpretation:** A significant p-value does not mean the behavior
  is *causal*; it only means the observed degree of synchrony is
  statistically unlikely to have occurred if the two participants were
  moving independently. For a pseudo-dyad null, it means more synchrony
  than with a stranger in the same setting.

## References

Kleinbub, J. R., & Ramseyer, F. T. (2020). rMEA: An R package to assess
nonverbal synchronization in motion energy analysis time-series.
*Psychotherapy Research*.
<https://doi.org/10.1080/10503307.2020.1844334>

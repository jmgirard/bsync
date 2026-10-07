# RB01: Segment-surrogate size failure (M015)

- **Date:** 2026-10-07
- **Output required:** write findings to `cairn/reviews/RR01-segment-surrogate-size.md`
- **Binding criteria:** not requested

You are performing an independent expert review. This brief is fully
self-contained. Do not assume any conversation context. Read only what this
brief directs you to read, answer the numbered questions, and write your
findings to the output path above using the same numbering.

## Background

`bsync` is an R package (C++ cores through Rcpp/RcppArmadillo) for
interpersonal synchrony in dyadic time series. Its windowed estimators are
windowed cross-correlation (`wcc()`), windowed dynamic time warping
(`wdtw()`), windowed Granger causality (`wgranger()`), and windowed phase
synchrony (`wphase()`). Each has a surrogate wrapper (`wcc_surrogate()` and
so on). A wrapper computes one aggregate statistic on the observed pair
`(x, y)` and the same statistic on `(x, y_surrogate)` for each column of a
user-supplied surrogate matrix. It reports the add-one Monte Carlo p-value
p = (b + 1) / (n + 1), where b counts surrogate statistics at least as
extreme as the observed one and n is the number of columns (Phipson & Smyth
2010, p. 6; `cairn/references/phipson2010.md`). A result is significant at
p <= .05. The default WCC statistic, `"mean_abs_z"`, is the mean of
|Fisher z| over every window and lag of the surface.

Milestone M015 adds `generate_surrogate_segment(y, segment_size,
n_surrogates)`. It follows the segment convention of the SUSY package
(Tschacher & Meier 2020; `cairn/references/tschacher2020.md`): it cuts `y`
into k = floor(length(y) / segment_size) segments from sample 1 and leaves
the tail (the last length(y) - k * segment_size samples) in place. As
planned and implemented, each surrogate column puts the k segments in an
order where no segment keeps its original position (a derangement). Orders
are distinct. With up to 8 segments, the generator lists all derangements
and draws distinct rows. With more, it draws permutations and rejects fixed
points and repeats. If `n_surrogates` is at least the number of
derangements, it returns all of them. The plan also adds
`surrogate_method = "segment"` to `synchrony_multiverse()` and
`autotune_wcc()`, with segment size = each grid cell's window size. In the
multiverse, the window increment is a fraction of the window
(`increment_pct`, default 0.1), so windows overlap.

The plan's size check was: independent pairs of AR(1) series
x_n = 0.7 x_{n-1} + eta_n, `wcc_surrogate()` with 19 segment surrogates,
window_size = window_increment = segment_size, lag_max >= 1, at least 6
segments, rejection rate at p <= .05 at most 0.09.

That check failed. Measured rates (1000 independent pairs per row, 19
surrogates, window 32, lag_max 4, `"mean_abs_z"`, nominal 0.05, Monte Carlo
SE about 0.007):

| Orders drawn | Segment size | Increment | Segments k | Rate |
|---|---|---|---|---|
| derangements | 32 | 32 | 8 | 0.090 |
| derangements | 32 | 32 | 5 | 0.118 |
| derangements | 32 | 32 | 20 | 0.073 |
| derangements | 32 | 3 | 8 | 0.104 |
| derangements | 40 | 40 | 8 | 0.063 |
| derangements | 40 | 3 | 8 | 0.073 |
| any non-identity | 32 | 32 | 8 | 0.068 |
| any non-identity | 32 | 32 | 5 | 0.049 |
| any non-identity | 32 | 3 | 8 | 0.066 |
| any non-identity | 40 | 40 | 8 | 0.049 |
| any non-identity | 40 | 40 | 5 | 0.043 |
| any non-identity | 40 | 3 | 8 | 0.060 |
| phase randomization | — | 32 | — | 0.059 (lag_max 1, n = 256) |
| phase randomization | — | 3 | — | 0.055 (n = 640) |

Also: derangements with segment = window = increment = 32 and lag_max 1
gave 0.074 (k = 8) and 0.083 (k = 20).

The implementing session proposed two causes:

1. Derangements are not a group that contains the identity. Under the null,
   the observed assignment (the identity) shares no position with any
   surrogate, but two surrogates share about one position on average, so
   the surrogate statistics are more alike than the observed one is to
   them. Drawing any non-identity permutation (fixed segments allowed)
   restores the group structure.
2. Lagged windows cross segment boundaries. In `wcc()`, window r of `x`
   covers samples `1 + lag_max + (r - 1) * increment` onward, and the `y`
   window at lag tau is shifted by tau. With segment = window, most lagged
   `y` windows straddle two segments. The observed `y` is continuous across
   that boundary and the surrogate is not, which lowers the within-window
   autocorrelation of surrogate `y` and so the spread of its correlations.
   With segment_size = window_size + 2 * lag_max and
   window_increment = segment_size, every lagged `y` window of window r
   lies inside segment r.

The proposed fix is both changes. Its rates are 0.043 to 0.049 in the
aligned design and 0.060 to 0.066 with overlapping windows. The plan's
Goal names the derangement rule, and the user chose window size as the
multiverse segment size, so the fix needs a re-plan. The user asked for an
independent review before that.

## Materials

Read these:

- `R/surrogate_generation.R` lines 340-480: `generate_surrogate_segment()`,
  `.n_segment_orders()`, `.segment_orders()` as implemented.
- `R/surface.R` lines 94-123: `build_surface_grid()`, which places windows
  (`i_vals <- 1L + tau_max + (grid$row - 1L) * w_inc`) and lags.
- `src/wcc.cpp` lines 1-110: the window and lag indexing (`x[j]` with
  `y[j + tau]`).
- `R/surrogate_analysis.R` lines 1-140 (`wcc_surrogate()`), and lines
  380-430 (`wgranger_surrogate()`, which builds a non-lagged grid with
  `lagged = FALSE`).
- `R/multiverse.R` lines 280-420: surrogate matrices built once per method
  and the per-cell loop (window, lag, and increment in samples).
- `cairn/references/tschacher2020.md` (SUSY segment convention: SUSY
  correlates segment pairs i != h and never builds a reordered series) and
  `cairn/references/phipson2010.md` (add-one p-value).
- `cairn/milestones/M015-segment-shuffle-surrogates.md`: Scope, acceptance
  criteria, and the work log (the 2026-10-07 T3 lines).
- `tests/testthat/test-surrogate-calibration.R`: the existing IAAFT size
  test, as a model of the test style.

To reproduce or extend the table, run this script from the repo root with
`Rscript <file> <window> <lag_max> <increment> <n> <segment> <mode>`, where
mode is `derange` (the current generator) or `any` (any non-identity
order). It takes about one minute per row. Running new rows is welcome.

```r
devtools::load_all(".", quiet = TRUE)
ar1 <- function(n, phi = 0.7) {
  as.numeric(stats::filter(stats::rnorm(n), phi, method = "recursive"))
}
# Any non-identity order, distinct, without replacement.
seg_any <- function(y, s, m) {
  n <- length(y); k <- n %/% s
  seen <- character(0); out <- matrix(0, n, 0)
  while (ncol(out) < min(m, factorial(k) - 1)) {
    p <- sample.int(k); key <- paste(p, collapse = ",")
    if (all(p == seq_len(k)) || key %in% seen) next
    seen <- c(seen, key)
    rows <- c(rep((p - 1) * s, each = s) + seq_len(s), seq_len(n - k * s) + k * s)
    out <- cbind(out, y[rows])
  }
  out
}
a <- commandArgs(TRUE); num <- as.numeric(a[1:5])
w <- num[1]; lag <- num[2]; inc <- num[3]; n <- num[4]; seg <- num[5]; mode <- a[6]
set.seed(777)
p <- vapply(seq_len(1000), function(i) {
  x <- ar1(n); y <- ar1(n)
  surr <- if (mode == "any") seg_any(y, seg, 19) else
    suppressMessages(generate_surrogate_segment(y, seg, 19))
  suppressWarnings(wcc_surrogate(x, y, surr, window_size = w, lag_max = lag,
    window_increment = inc)$p_value)
}, numeric(1))
cat(sprintf("%s w=%d lag=%d inc=%d n=%d seg=%d k=%d rate=%.4f\n",
  mode, w, lag, inc, n, seg, n %/% seg, mean(p <= .05)))
```

## Questions

1. Is cause 1 correct? State the exact condition under which the add-one
   p-value of a segment-permutation test is valid (size at most alpha)
   under the null of independent `x` and `y`, and whether sampling
   derangements breaks it. Is sampling distinct non-identity permutations
   without replacement from all k! - 1, with p = (b + 1) / (m + 1), exact
   under that condition (phipson2010)? Does this still hold when m is
   below k! - 1, and when all k! - 1 are used?
2. Is cause 2 correct, given the window and lag indexing in `surface.R` and
   `wcc.cpp`? Is segment_size = window_size + 2 * lag_max with
   window_increment = segment_size the right alignment condition for
   `wcc()`? What is the condition for `wdtw()` and `wphase()` (lagged grids)
   and for `wgranger()` (non-lagged grid, AR order `ar_order`)?
3. With both fixes, what residual departures from exactness remain (for
   example, correlation between adjacent segments of an autocorrelated `y`,
   the tail kept in place, unused samples between windows)? Which null
   hypothesis does the fixed test then test, stated so a user can read it?
4. The multiverse uses overlapping windows (increment 10% of the window by
   default). Measured size with both fixes was about 0.060 to 0.066. Which
   option do you recommend for bsync: (a) for `"segment"` cells, force
   window_increment = segment_size = window_size + 2 * lag_max; (b) keep
   the cell's increment and document the measured excess; (c) do not offer
   `"segment"` in the multiverse and `autotune_wcc()` at all; (d) another
   design? Weigh validity against the user's stated wish to offer the
   method in both functions.
5. Do you recommend keeping the derangement form available at all, for
   example as an option for users who want SUSY's "non-matching segments"
   semantics? If bsync drops it, what do the docs need to say about SUSY's
   own surrogate comparison?
6. Propose a size test design for the fixed generator: process, series
   length, k values (including small k, where k! - 1 can be below 19),
   number of pairs, number of surrogates, and rejection-rate bounds with
   the Monte Carlo error that justifies them. Include a control that shows
   the test can fail (a planted defect such as derangement-only sampling).

## Constraints

- D-001: p-values use the add-one form (b + 1) / (n + 1). Fixed.
- D-002: significance is p <= .05. Fixed.
- Invariant 2 (matched null): the observed and surrogate statistics must be
  computed identically. Fixed.
- Each surrogate wrapper needs a surrogate matrix with exactly `length(y)`
  rows, so a generator that drops samples must say how it fills them.
- The plan gate chose a reordered full-length series over SUSY-exact segment
  pairing (correlating segment pairs i != h), because the reordered series
  works with all four wrappers and the multiverse. Flag explicitly if you
  think that choice is unsound rather than working around it.
- Do not edit any repo file except the output RR.

## Output format

In `RR01-segment-surrogate-size.md`: answer each question by number with
your reasoning and evidence; list any additional findings separately under
"Beyond the brief"; end with concrete recommendations, each marked apply /
consider / reject-with-reason. Your report is advisory: emit a
`## Binding criteria` section ONLY if this brief's header slot says
`requested`.

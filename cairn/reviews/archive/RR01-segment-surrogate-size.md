# RR01: Segment-surrogate size failure (M015)

- **Date:** 2026-10-07
- **Brief:** `cairn/reviews/RB01-segment-surrogate-size.md`
- **Binding criteria:** not requested (none emitted)

## Materials and new evidence

Read: `R/surrogate_generation.R` 340-480. `R/surface.R` 94-160. The index
lines of `src/wcc.cpp`, `src/wdtw.cpp` and `src/wphase.cpp`.
`src/wgranger.cpp` 1-100. All four wrappers in `R/surrogate_analysis.R`.
`R/multiverse.R` 280-660. `R/autotune.R` 134-260. Both reference notes.
The SUSY 0.1.0 shelf copy `cairn/references/sources/susy-0.1.0-susy.R`,
lines 95-165 and 305-330, for Q5. The M015 milestone file.
`tests/testthat/test-surrogate-calibration.R`.

Two simulation sets were run. The scripts are in the session scratchpad
and are not committed.

Set A is an idealized model with no time-series structure. It has `k`
independent white-noise x-windows and `k` independent white-noise
y-segments of length 32. The statistic is the mean of |Fisher z| over the
paired windows. There are 19 surrogates and 20000 reps per row. Its SE is
0.0015.

Set B is the brief's script with four additions: the wrapper, the AR
coefficient, a `rotate` mode, and a seed. Each row has 1000 independent
AR(1) pairs unless stated, 19 surrogates, and the add-one p-value at
p <= .05. Its SE is 0.007. "Aligned" below means
`segment_size = window_increment = window_size + 2 * lag_max`. The full
Set B table is in the appendix.

## 1. Cause 1 (derangements) is correct

**Exact condition.** Write the segments of `y` as S_1..S_k. Under the
null, `x` and `y` are independent. The permutation test has size at most
alpha under two conditions.

(i) *Invariance.* Permuting the segments of `y` leaves the joint law of
`(x, y)` unchanged: `(x, y_g)` has the law of `(x, y)` for every `g` in a
group `G` that acts on segment positions. Because `x` is independent of
`y`, this is a statement about `y` alone. Its segments must be
exchangeable (then `G` is the full symmetric group `S_k`), or at least
invariant under `G`.

(ii) *Symmetric sampling.* The observed arrangement must not be
distinguishable from the surrogates by the shape of the sampled set. A
sufficient rule: draw `m` distinct elements uniformly without replacement
from `G \ {id}`, where `G` is a group that contains the identity.

Proof sketch of (ii). Let `Sigma` be uniform on `G` and independent of the
data. By (i), the observed data has the law of `(x, y_Sigma)`. The
surrogates are then `sigma_j Sigma`. The test evaluates the statistic on
the set `U = {Sigma, sigma_1 Sigma, ..., sigma_m Sigma}`, which lies in
the coset `G Sigma`. For a given `U`, `P(Sigma = g | U)` is proportional
to `P(Sigma = g)` times the chance that `{u g^-1 : u in U, u != g}` was
the sampled set. For every `g` in `U`, that set is an `m`-subset of
`G \ {id}`, and every such subset is equally likely. So `Sigma` is uniform
within `U`. Its rank among the `m + 1` statistics is uniform. Ties count
against rejection under the `>=` rule. So
`P(p <= alpha) = floor(alpha (m + 1)) / (m + 1) <= alpha`. This is the
p. 6 result of Phipson and Smyth, in which `m` counts the non-identity
draws and the `+ 1` is the observed arrangement.

**Derangements break (ii).** The set `D_k` of derangements is not a
group. It lacks the identity, and it is not closed under composition.
Take `g = sigma_1 Sigma` in `U`. The relabeled set contains
`sigma_j sigma_1^-1`. That element is a derangement only for a pair
`sigma_j`, `sigma_1` that agrees at no position. Two independent uniform derangements agree
at `k / (k - 1)` positions on average, about 1. Each surrogate agrees with
the observed arrangement at 0 positions. So the observed arrangement is
identifiable as the odd one out. Every surrogate shares no term of the
sum with it. Surrogates share about one of `k` terms with each other. The
surrogate statistics are therefore more alike than the observed statistic
is to them. The observed one is extreme more often than `1 / (m + 1)` of
the time. The effect shrinks like `1 / k`.

Under "any non-identity" the shape is symmetric. A uniform permutation has
one fixed point on average, and the proof above does not depend on the
fixed-point counts.

**Evidence (Set A).**

| k | derangements | any non-identity | rotations (k - 1 shifts) |
|---|---|---|---|
| 4 | 0.000 (all 9 orders, p_min = 0.10) | 0.048 | none |
| 5 | 0.113 | 0.052 | none |
| 8 | 0.080 | 0.047 | none |
| 12 | 0.069 | 0.049 | none |
| 20 | 0.060 | 0.053 | 0.051 |

The derangement inflation appears with no autocorrelation, no lags and no
windows. It is a property of the sampling set alone. The `k = 5, 8, 20`
rows match the brief's aligned derangement rows (0.118, 0.090, 0.073, of
which the last two also include cause 2). They also match Set B's aligned
derangement rows: `k = 5` gave 0.089 to 0.102 over four seeds, `k = 8`
gave 0.080, `k = 12` gave 0.060.

**Is "any non-identity, distinct, add-one" exact?** Yes, under (i), for
every `m` from 1 to `k! - 1`.

- For `m < k! - 1` it is exact by the argument above (Phipson and Smyth
  section 4, p. 6).
- For `m = k! - 1` (all orders), `U = G` and `p = (b + 1) / k!`. This is
  the exhaustive permutation p-value (section 5, p. 7). Exact.
- Draws *with* replacement from `G` make `(b + 1) / (m + 1)` conservative
  (section 6.1, p. 7). The plan's without-replacement choice is right.

Two consequences for small `k`. With `k = 3` there are 5 non-identity
orders, so `p_min = 1/6` and the test can never reject at .05. With
`k = 4` there are 23, so `p_min = 1/24` and 19 of them can be drawn.
Rejection at .05 needs `k! >= 20`, that is `k >= 4`. Set B at `k = 4`
(19 of 23 orders) gave 0.058.

## 2. Cause 2 (lagged windows cross segment boundaries) is correct

**Indexing.** `build_surface_grid()` places x-window `r` (1-based) at
`a_r = 1 + L + (r - 1) * inc`. The window covers `x[a_r .. a_r + w - 1]`.
The C++ cores pair `x[j]` with `y[j + tau]`. So the y-window at lag `tau`
covers `y[a_r + tau .. a_r + tau + w - 1]`. Over `tau = -L..L`, the
y-windows of window `r` together span `y[1 + (r - 1) inc .. (r - 1) inc +
w + 2L]`. That is a block of `w + 2L` samples that starts at
`1 + (r - 1) inc`.

Now set segment = window = increment = `s = w` with `L >= 1`. Segment `r`
is `[(r - 1) s + 1, r s]`. The y-window at `tau = -L` is that segment
exactly. Every other lag straddles segments `r` and `r + 1`. So `2L` of
the `2L + 1` lags straddle. In the observed series that stretch is one
continuous AR(1) run. In a surrogate it has a jump between two unrelated
segments. The jump lowers the within-window serial dependence of surrogate
`y`. That lowers the spread of the window correlation with an
autocorrelated `x`, which lowers mean |z|. The observed value then sits
high in the null. Set B isolates this with the fixed sampling rule:
`w = 16, L = 2, s = inc = 16` gave 0.090, and `s = inc = 20 = w + 2L`
gave 0.050.

**Alignment for `wcc()`.** The block above lies inside segment `r` under
two conditions: `inc = s` and `s >= w + 2L`. The case `s = w + 2L` is
tight. Any larger `s` with `inc = s` also works and leaves `s - w - 2L`
unused samples at the segment ends (a buffer, see Q3). `lag_increment`
does not matter. So the brief's condition is right. It is the tight case
of `segment_size >= window_size + 2 * lag_max` with
`window_increment = segment_size`.

**`wdtw()` and `wphase()`.** Both use the same lagged grid. Both pair
`x[i + w]` with `y[i + tau + w]` (`wdtw.cpp` 29, 45, 90, and `wphase.cpp`
67, 82). So the same condition holds at the window level.
Two qualifications apply.

- `wdtw()` with `scale_method = "global"` z-scores each surrogate column
  by its own mean and SD. A segment permutation leaves both unchanged. So
  the surrogate window contents are exactly the segment contents of the
  scaled `y`, and the alignment condition is sufficient. Set B (`w = 16,
  L = 2, s = inc = 20, k = 5`, 500 pairs) gave 0.044. The
  `fast_method = TRUE` path uses `tau = 0` over the same windows, a
  sub-block, so it is covered.
- `wphase()` applies `gsignal::hilbert()` to the whole series before
  windowing. The analytic signal of a segment-shuffled series is not the
  segment-shuffle of the analytic signal. The `1 / (pi t)` kernel is
  global, and the jumps at segment joins ring into the nearby samples. So
  window alignment does not make the phase test exact. The surrogate
  phases inside segment `sigma(r)` differ from the observed phases inside
  that segment near its ends. Set B (`w = 32, L = 4, s = inc = 40,
  k = 8`) gave 0.059 on AR(1) input. That is 1.3 SE above nominal and
  settles nothing. An exact form permutes the phase series `phi_y` by segments
  instead of `y`. The current generator and wrapper split cannot express
  that. Treat segment surrogates for `wphase()` as approximate until they
  are measured on an oscillatory process, and say so in the docs.

**`wgranger()`.** The non-lagged grid places window `r` at
`1 + (r - 1) inc`. The core uses only `x[i .. i + w_max]` and
`y[i .. i + w_max]`. The regression rows are `t = i + p + row` with lags
`t - 1 .. t - p`, all inside the window (`wgranger.cpp` 58-74). The
condition is `segment_size = window_size` with
`window_increment = window_size`. `ar_order` does not enter it. Set B
(`w = inc = s = 32, k = 8`, 5000 pairs over three seeds) gave 0.062,
0.051 and 0.044, pooled 0.050. The derangement control at the same design
gave 0.103. Straddling matters little for Granger: `s = 40, inc = 32`
gave 0.057 and the `inc = 3` overlap gave 0.044. That fits an F statistic
that depends on lag-1 structure, which one jump barely disturbs.

**Window count.** In the aligned design with `n = k s` exactly, the grid
yields `k - 1` windows, not `k` (see Beyond the brief, B1). Segment `k`
then appears only as a donor in surrogates. Exactness is unaffected,
because the proof in Q1 holds for any fixed statistic. It costs one
window, and the docs need to state it.

## 3. Residual departures with both fixes, and the null tested

With "any non-identity" sampling and the aligned design, the add-one
p-value is exact under (i), the segment exchangeability of `y`. What
remains is the failure of (i).

1. **Serial dependence across segment joins.** The segments of an
   autocorrelated `y` are not exchangeable. The first samples of segment
   `r + 1` depend on the last samples of segment `r`. In the observed
   arrangement, adjacent x-windows and adjacent y-windows are both
   time-adjacent. So neighbouring window statistics covary by about
   `sum_d d rho_x(d) rho_y(d) / s^2`. In a surrogate the y-neighbours are
   unrelated, and that covariance is near zero. The observed mean |z| is
   therefore a little more variable than the surrogates' values, which
   inflates both tails. Relative to the variance of a window statistic,
   the covariance is about `rho^2 / (s (1 - rho^4))`: 0.016 for
   `rho = 0.7, s = 40`, and 0.12 for `rho = 0.95`. Set B at `rho = 0.95`
   (`s = 20`, `k = 5` and `k = 12`) gave 0.047 and 0.054. So the effect is
   not visible at 1000 pairs even for a strongly autocorrelated AR(1). It
   grows for long-memory or slowly drifting `y`. The aligned design
   already reserves `2L` samples per segment outside the x-window. A
   `segment_size` above `w + 2L` adds a buffer that weakens this
   dependence further.
2. **Tail kept in place.** In the aligned design no window touches the
   tail. For tail length `t`, `n_r = k - 1 + floor((t + 1) / s)`, so the
   last window ends at or before sample `k s`. The tail is harmless there.
   With overlapping windows, the tail is in place in the observed series
   and in every surrogate alike, so it is symmetric. The one asymmetry is
   the join between the last segment and the tail. Windows that straddle
   it are continuous in the observed series and not in surrogates, the
   same mechanism as cause 2.
3. **Unused samples.** Segment `k` (when `n_r = k - 1`), the tail, and
   any buffer samples contribute to surrogates as donors but never to the
   observed statistic. Exactness holds. Power is a little lower than the
   full series gives.
4. **`x` is not permuted.** Only the segments of `y` need to be
   exchangeable. The x-windows are fixed and can be nonstationary without
   harm. For `wcc()` the Pearson correlation is invariant to each window's
   mean and scale. So level changes between segments of `y` do not break
   (i) for `wcc()`. They do break it for `wdtw()` (global scaling) and for
   the analytic signal in `wphase()`.
5. **Phase extraction** is the one wrapper-level departure that alignment
   does not remove (Q2).

**The null, stated for a user.** "Cut `y` into blocks of `segment_size`
samples. Each window of `x` has one block of `y` that occurred alongside
it. The test asks whether that block resembles the window more than the
other blocks of `y` do. The null hypothesis says it does not. Under the
null, the time at which each block of `y` occurred carries no information
about `x`. Any arrangement of the blocks is then as good as the real one.
Everything inside a block is kept: its autocorrelation, its distribution,
and any trend within it. Only the co-occurrence of `x` and `y` at the
block scale is destroyed. The test is exact under two conditions. First,
the blocks of `y` must be exchangeable. This is nearly true for a block
that is long relative to the memory of `y`. Second, every window and every
lagged window must lie inside one block. This holds for
`segment_size >= window_size + 2 * lag_max` with
`window_increment = segment_size`."

## 4. The multiverse: option (a), with the alignment stated loudly

Measured excess with overlapping windows (10% increment): 0.060 to 0.066
in the brief. Set B gave 0.067 (`w = 16, L = 2, s = 20, inc = 2,
k = 12`). Three rows, each 1.5 to 2.4 SE above nominal, all above. That
is a real excess of about +0.015 at nominal .05. Phase randomization in
the same design gave 0.055 to 0.059. So the segment method with overlap is
about one SE worse than the package's default method.

- **(a) Force `window_increment = segment_size = window_size + 2 lag_max`
  for `"segment"` cells.** This restores the exact design. The costs: the
  cell's `increment_pct` is overridden. The observed statistic becomes a
  mean over `k - 1` non-overlapping windows rather than over many
  overlapping ones (the same estimand, a noisier estimate). And
  `n_windows` differs from the other methods' cells at the same window
  and lag. This is how SUSY itself works. Its segments are its windows,
  they do not overlap, and lags are taken inside the segment (`ccf()` on
  the two segments). The method's own lineage is non-overlapping and
  aligned.
- **(b) Keep the cell's increment and document the excess.** This honours
  the increment axis but ships a test with a documented size of about
  0.065 at .05. "Loud defaults" is met by the documentation. "Safe
  defaults" is not. The existing IAAFT calibration bound of 0.08 admits
  it, which says more about that bound than about this method.
- **(c) Do not offer `"segment"` in the multiverse.** This rejects the
  user's stated wish and gains nothing over (a).
- **(d) Other designs.** Restricting the statistic to the windows that
  lie inside one segment collapses to (a) at `s = w + 2L`. Per-window
  pairing (SUSY-exact) was rejected at the plan gate and does not fit the
  shared engine. A segment-permuted phase series for `wphase()` is the
  only (d) worth a candidate row.

**Recommendation: (a).** Implement it as follows. For each `"segment"`
cell set `s = w_samp + 2 * l_samp` (`s = w_samp` for Granger) and
`inc_samp = s`. Generate one surrogate matrix per distinct `s`. Record the
realized `window_increment` in the grid. Skip cells with
`s > floor(n / 2)` through the existing skip path, with a message that
names the window and lag. Emit one `cli_inform()` per call that says
segment cells use aligned non-overlapping windows and why. In the
`synchrony_multiverse()` and `autotune_wcc()` docs, state that
`increment_pct` does not apply to segment cells and that their
`n_windows` is about `floor(n / s) - 1`. The user can prefer
comparability of the observed statistic across methods over exact size.
In that case (b) is defensible only with the measured excess printed in
the roxygen and the vignette. I do not recommend it.

Note for `autotune_wcc()`. Under (a), the number of orders `k! - 1` and
the number of windows both fall as the window grows. The default
`n_surrogates = 100` is reached only for `k >= 5` (`5! - 1 = 119`). With
`k = 4` the generator returns all 23 orders. `warn_few_surrogates()`
tests the requested count, not the delivered one. The cell's `p_min` needs
the delivered count.

## 5. Keep derangements? No. Rotations are the exact substitute

Derangement sampling is anti-conservative by construction (Q1). It has no
compensating property for a permutation p-value. Its one legitimate use
is a null-mean and SD reference in which no surrogate carries a true-pair
segment. That raises power by about `1 / k` relative to uniform
permutations. The same property is available exactly through a
**fixed-point-free group**: cyclic rotations of the segments,
`sigma_j(r) = r + j mod k` for `j = 1..k - 1`. Every non-identity rotation
moves every segment. The rotations with the identity form a group, so
sampling distinct rotations is exact by Q1. Set A at `k = 20` gave 0.051.
Set B through `wcc_surrogate()` at `k = 20` gave 0.055. The cost: only
`k - 1` surrogates exist, so `p_min = 1 / k`, and 19 surrogates need
`k >= 20`. Offer it as an option only on request. The candidate-row cost
is small, and the M015 cost is not zero. Drop derangements from the
package.

**What the docs must say about SUSY once derangements go.** From the
shelf copy (lines 95-165 and 305-330): SUSY never builds a reordered
series. It cross-correlates every segment pair `(i, h)` of the two series
with `ccf()` up to `maxlag`. Pairs with `i == h` are real, pairs with
`i != h` are pseudo. It averages |Fisher z| per lag over the real pairs
and over all `k (k - 1)` pseudo pairs. It reports
`ES = (mean_lags realZ - mean_lags pseudoZ) / sd_lags(pseudoZ)` and the
percentage of lags at which real exceeds pseudo. It computes no
permutation p-value. So "no segment in place" describes SUSY's *pair
set*, not any surrogate series, and the derangement rule was bsync's
translation of it. The docs need four statements. First, bsync's segment
surrogates are reordered full-length series drawn uniformly from all
non-identity orders. A surrogate can therefore leave about one segment of
`k` in place. Second, this is required for the add-one p-value to be exact, and
it is why bsync does not follow SUSY's non-matching rule. Third, bsync
does not reproduce SUSY's ES or its numbers. Fourth, SUSY's lag rule
`segment >= 2 * maxlag` has the bsync analogue
`segment_size >= window_size + 2 * lag_max`.

## 6. Size test design for the fixed generator

Process: independent AR(1) pairs with `phi = 0.7`, as in the IAAFT test.
Wrapper: `wcc_surrogate()`, `"mean_abs_z"`, 19 surrogates. Under
exactness `P(p <= .05) = 1/20 = 0.05` exactly, and the binomial SE at
1000 pairs is 0.0069.

Design: `window_size = 16`, `lag_max = 2`, `segment_size =
window_increment = 20`. Two positive runs cover both code paths.

- Run 1: `k = 5`, `n = 100` (enumeration path, `5! - 1 = 119` orders).
- Run 2: `k = 12`, `n = 240` (rejection-sampling path).

Each run has 1000 pairs and takes about 22 s on this machine. Assert
`rate <= 0.07` and `rate >= 0.03`. The upper bound is 2.9 SE above 0.05,
so a false failure under exactness has a chance of about 0.002 per run.
The lower bound guards against a collapsed null. Over eleven seeds of the
aligned fixed design in Set B, the rates ran from 0.044 to 0.058.

Control (planted defect): Run 1's design with derangement-only sampling,
1000 pairs. Assert `rate > 0.07`. Four seeds gave 0.089, 0.098, 0.099 and
0.102. At a true rate of 0.095, the chance of falling to 0.07 or below is
about 0.004. The control must use small `k`. At `k = 12` the derangement
rate is 0.060, and the bound does not catch it. Build the planted defect
in the test file in one of two ways. Filter the fixed generator's output
to columns with no segment in place, or write a small reference
generator. The IAAFT
test's "values form" control is the model.

Small `k`: add a unit test that `k = 3` returns all 5 non-identity orders
and that the wrapper's `p_min` is `1/6`. Add one that `k = 4` with 19
requested returns 19 distinct non-identity orders. Do not run a size
check at `k = 3`, because its rate is 0 by construction. Replace AC3's
"second run with k = 3" with the `k = 4` or `k = 5` run above.

Optional extension rows, not for the test suite: the tail case
(`n = 113, k = 5`: 0.047), `rho = 0.95` (0.047 and 0.054), and Granger
aligned (`s = inc = w = 32`, pooled 0.050).

## Beyond the brief

- **B1. `build_surface_grid()` drops one valid window when
  `window_increment > 1`.** Let `q` be `n` minus `w_max` minus `2L`. The
  code sets `n_r` to the floor of `q / inc`. The correct count is the
  floor of `(q - 1) / inc`, plus 1. When `inc` divides `q`, the two agree.
  Otherwise the code is one window short. Checked:
  `build_surface_grid(5, 3, 5, 1)` aborts with "Need at least 6 samples",
  although `x[2..4]` with `y[1..5]` fits in 5. And
  `build_surface_grid(k * 40, 32, 40, 4)$n_r` is `k - 1`. The non-lagged
  branch has the same form, and the multiverse's `min_n` guard copies it.
  This is pre-existing and outside M015. `inc = 1` is unaffected, so no
  shipped default is wrong. It costs the aligned segment design one
  segment. A candidate row, not a hotfix.
- **B2. Significance is unreachable below `k = 4`.** `k! - 1 < 19` for
  `k <= 3`. The generator or the wrapper needs to say so, as
  `print_surrogate_count_notes()` already does for few surrogates. The
  multiverse needs `p_min` from the delivered column count.
- **B3. `.n_segment_orders()` and the above-8 abort.** Under the fix the
  order count is `k! - 1`, so `9! - 1 = 362879`. The abort "must be below
  the possible orders" is then reachable only for absurd `n_surrogates`.
  Keep the enumeration path for `k <= 8` and the all-orders return. Above
  8, reject the identity and repeats.
- **B4. Phase randomization at 0.055 to 0.059.** Two rows in the brief,
  both above nominal by about 1 SE. Not alarming. But the package has no
  recorded size for its default method on a plain AR(1) pair. One
  2000-pair row settles whether this is noise.
- **B5. The IAAFT calibration bound (0.08 at 400 pairs).** It is 2.7 SE
  above nominal and admits the overlapping-window segment design the brief
  measured. Tighten it to 0.07 at 1000 pairs at the next edit of that
  file.

## Recommendations

1. **Apply.** Replace derangement sampling with uniform distinct
   non-identity orders in `generate_surrogate_segment()`. Keep the
   enumeration path for `k <= 8` and the rejection path above. If `n_surrogates`
   reaches `k! - 1`, return all orders. Update the
   Goal, Scope and the AC1/AC4 wording ("no segment in place" goes).
2. **Apply.** Document the alignment condition. For `wcc()`, `wdtw()`
   and `wphase()`: `segment_size >= window_size + 2 * lag_max` with
   `window_increment = segment_size`. For `wgranger()`:
   `segment_size = window_size = window_increment`. Give the straddling
   mechanism in one sentence and the Q3 null statement.
3. **Apply.** Multiverse option (a). Per segment cell `s = w + 2 L`
   (`s = w` for Granger) and `inc = s`. One matrix per distinct `s`. The
   realized increment in the grid. One informational message per call.
   If `s > floor(n / 2)`, skip the cell.
4. **Apply.** The Q6 test. Two positive runs (`k = 5`, `k = 12`), bounds
   [0.03, 0.07] at 1000 pairs, a derangement control at `k = 5` asserted
   above 0.07, and the small-`k` unit tests. Drop the `k = 3` size run.
5. **Apply.** Docs on SUSY as in Q5: pair set versus reordered series, no
   SUSY numbers reproduced, the lag-rule analogue.
6. **Consider.** Segment rotations (`k - 1` exact, fixed-point-free
   surrogates) as a candidate row for users who want every surrogate to
   move every segment.
7. **Consider.** Flag `wphase()` segment surrogates as approximate. Add a
   candidate row to measure them on an oscillatory process, or to permute
   `phi_y` by segments inside `wphase_surrogate()`.
8. **Consider.** B1 (grid window count) and B2 (unreachable significance
   at `k <= 3`) as candidate rows. B4 and B5 as a tidy-up at the next edit
   of `test-surrogate-calibration.R`.
9. **Reject.** Option (b) for the multiverse. It ships a known +0.015
   size excess in a default, and the only gain is cross-method
   comparability of `n_windows`.
10. **Reject.** Option (c). The aligned design makes the method sound in
    the multiverse, and withholding it loses the user's stated goal.
11. **Reject.** Keeping derangements as a p-value option, even flagged.
    It is anti-conservative at every `k` and never exact.

The plan gate's choice of a reordered full-length series over SUSY-exact
pairing is sound. With aligned windows and uniform non-identity sampling
it is an exact permutation test, and the shared engine can run it for
every wrapper. Pairing cannot do that.

## Appendix: Set B rows

Format: wrapper, sampling, w, lag, inc, n, segment, k, phi, pairs, rate
at p <= .05 (rate at p <= .10).

```
wcc      any     w=16 lag=2 inc=20 n=80  seg=20 k=4  phi=0.70 pairs=1000 rate=0.058 (0.113)
wcc      any     w=16 lag=2 inc=20 n=100 seg=20 k=5  phi=0.70 pairs=1000 rate=0.046 (0.088)
wcc      any     w=16 lag=2 inc=20 n=100 seg=20 k=5  phi=0.70 pairs=1000 rate=0.044 (0.095)
wcc      any     w=16 lag=2 inc=20 n=100 seg=20 k=5  phi=0.70 pairs=1000 rate=0.056 (0.115)
wcc      any     w=16 lag=2 inc=20 n=100 seg=20 k=5  phi=0.70 pairs=1000 rate=0.055 (0.108)
wcc      any     w=16 lag=2 inc=20 n=113 seg=20 k=5  phi=0.70 pairs=1000 rate=0.047 (0.089)  tail 13
wcc      any     w=16 lag=2 inc=20 n=240 seg=20 k=12 phi=0.70 pairs=1000 rate=0.050 (0.096)
wcc      any     w=16 lag=2 inc=20 n=240 seg=20 k=12 phi=0.70 pairs=1000 rate=0.053 (0.108)
wcc      any     w=16 lag=2 inc=20 n=240 seg=20 k=12 phi=0.70 pairs=1000 rate=0.048 (0.093)
wcc      any     w=16 lag=2 inc=20 n=240 seg=20 k=12 phi=0.70 pairs=1000 rate=0.052 (0.100)
wcc      any     w=16 lag=2 inc=20 n=400 seg=20 k=20 phi=0.70 pairs=1000 rate=0.048 (0.083)
wcc      any     w=16 lag=2 inc=20 n=100 seg=20 k=5  phi=0.95 pairs=1000 rate=0.047 (0.093)
wcc      any     w=16 lag=2 inc=20 n=240 seg=20 k=12 phi=0.95 pairs=1000 rate=0.054 (0.095)
wcc      any     w=16 lag=2 inc=16 n=240 seg=16 k=15 phi=0.70 pairs=1000 rate=0.090 (0.151)  straddle only
wcc      any     w=16 lag=2 inc=2  n=240 seg=20 k=12 phi=0.70 pairs=1000 rate=0.067 (0.133)  10% overlap
wcc      rotate  w=16 lag=2 inc=20 n=400 seg=20 k=20 phi=0.70 pairs=1000 rate=0.055 (0.104)
wcc      derange w=16 lag=2 inc=20 n=100 seg=20 k=5  phi=0.70 pairs=1000 rate=0.102 (0.160)  control
wcc      derange w=16 lag=2 inc=20 n=100 seg=20 k=5  phi=0.70 pairs=1000 rate=0.099 (0.149)  control
wcc      derange w=16 lag=2 inc=20 n=100 seg=20 k=5  phi=0.70 pairs=1000 rate=0.089 (0.139)  control
wcc      derange w=16 lag=2 inc=20 n=240 seg=20 k=12 phi=0.70 pairs=1000 rate=0.060 (0.120)
wcc      derange w=32 lag=4 inc=40 n=200 seg=40 k=5  phi=0.70 pairs=1000 rate=0.098 (0.154)
wcc      derange w=32 lag=4 inc=40 n=320 seg=40 k=8  phi=0.70 pairs=1000 rate=0.080 (0.146)
wdtw     any     w=16 lag=2 inc=20 n=100 seg=20 k=5  phi=0.70 pairs=500  rate=0.044 (0.094)
wphase   any     w=32 lag=4 inc=40 n=320 seg=40 k=8  phi=0.70 pairs=1000 rate=0.059 (0.106)
wgranger any     w=32       inc=32 n=256 seg=32 k=8  phi=0.70 pairs=1000 rate=0.062 (0.111)
wgranger any     w=32       inc=32 n=256 seg=32 k=8  phi=0.70 pairs=2000 rate=0.051 (0.106)
wgranger any     w=32       inc=32 n=256 seg=32 k=8  phi=0.70 pairs=2000 rate=0.044 (0.097)
wgranger any     w=32       inc=32 n=240 seg=40 k=6  phi=0.70 pairs=1000 rate=0.057 (0.100)  misaligned
wgranger any     w=32       inc=3  n=256 seg=32 k=8  phi=0.70 pairs=1000 rate=0.044 (0.100)  overlap
wgranger derange w=32       inc=32 n=256 seg=32 k=8  phi=0.70 pairs=1000 rate=0.103 (0.160)  control
```

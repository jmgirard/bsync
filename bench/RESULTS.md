# bench/RESULTS.md — M2 WCC prefix-sum speedup

Measured on Apple clang 21 / macOS, R 4.6.1, 5 iterations each.
Hardware: Apple Silicon (arm64). `bench::mark()` median wall time.

## Before (O(n_c × n_r × w_max) loop, M1 baseline)

| Config           |     n | window_size | lag_max | Median (ms) |
|------------------|------:|------------:|--------:|------------:|
| sim_dyad_narrow  |  2400 |          96 |      40 |         270 |
| sim_dyad_wide    |  2400 |         240 |      80 |        1142 |
| large_narrow     | 10000 |         100 |      20 |         616 |
| large_wide       | 10000 |         500 |     100 |       13697 |

## After (O(n_c × n) prefix-sum, M2 rewrite)

| Config           |     n | window_size | lag_max | Median (ms) | Speedup |
|------------------|------:|------------:|--------:|------------:|--------:|
| sim_dyad_narrow  |  2400 |          96 |      40 |          50 |   **5.4×** |
| sim_dyad_wide    |  2400 |         240 |      80 |          93 |  **12.3×** |
| large_narrow     | 10000 |         100 |      20 |         113 |   **5.5×** |
| large_wide       | 10000 |         500 |     100 |         549 |  **24.9×** |

## Notes

Speedup scales with `window_size × lag_max` as expected: the prefix-sum
eliminates the per-pair inner loop (previously O(w_max) per (i, τ) cell),
replacing it with O(1) window evaluation after O(n) per-τ preprocessing.

The large_wide config (n=10,000, ws=500, lag=100) benefits most — 25× —
because it has the largest w_max (499) relative to n. The sim_dyad configs
show 5–12× because w_max is a smaller fraction of n=2400.

Numerical correctness confirmed by the `stats::cor` oracle at tolerance 1e-9
on sim_dyad (clean, NA na.rm=TRUE, NA na.rm=FALSE). Results are identical to
the prior implementation within floating-point precision.

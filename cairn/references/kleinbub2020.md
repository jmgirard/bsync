# kleinbub2020 — the rMEA shuffle() pseudo-dyad convention

**Provenance.** Ingested 2026-10-05 by M009 from the rMEA 1.2.2 CRAN source, read through the GitHub CRAN mirror (`gh api repos/cran/rMEA/contents/...`, mirror commit `22e8f87a35c3a7b271e327538f5931a5968d80c6`). Files read: `R/rMEA_rand.R` (`shuffle()`), `R/rMEA_util.R` (`unequalCbind()`), `DESCRIPTION`, and `inst/CITATION`. The article itself was not read.
Pagination: — (source code. Anchors are file and function names.)
Extraction: verified 2026-10-05 against the rMEA 1.2.2 source files named above, read directly — observed 2026-10-05.

**Citation.** Kleinbub, J. R., & Ramseyer, F. T. (2020). rMEA: An R package to
assess nonverbal synchronization in motion energy analysis time-series.
*Psychotherapy Research*, 1–14. DOI: 10.1080/10503307.2020.1844334. This is the
form that rMEA's `inst/CITATION` prints. The printed issue has a later year and
volume, which this note did not check.

**Role.** This note records the published pseudo-dyad convention that
`generate_surrogate_pseudo()` and `generate_pseudo_dyads()` follow, and the
places where bsync departs from it. The roxygen and the surrogate-testing
vignette cite the article.

## Extracted values

- Purpose. `shuffle()` "recombines the s1 and s2 motion energy time-series
  between all MEA objects in the supplied list", for comparing "genuine
  synchrony of real data with pseudosynchrony", `R/rMEA_rand.R`, roxygen.
- Pairing. The function pools all s1 series, then all s2 series, takes every
  unordered pair with `utils::combn()`, and removes each real pairing (s1 of
  dyad i with s2 of dyad i), `R/rMEA_rand.R`, `shuffle()` body.
- Roles. `keepRoles = FALSE` is the default. With `TRUE`, it removes s1–s1 and
  s2–s2 pairs, so every pseudo-dyad is an s1 of one dyad with an s2 of
  another. Pair counts that follow: N(N-1) with roles kept, and
  C(2N, 2) - N = 2N(N-1) with roles mixed.
- Sampling. `size = "max"` returns every pair. An integer `size` draws that
  many pairs with `sample()`, without replacement, `R/rMEA_rand.R`.
- Length. Each pair goes through `unequalCbind(..., keep = FALSE)`, which
  keeps rows `1:minlen` of each series (a start-aligned crop to the shorter
  series), then `stats::na.omit()`, `R/rMEA_util.R` and `R/rMEA_rand.R`.
- Within-subject segment shuffling (`shuffle_segments()`) is disabled in
  1.2.2. The roxygen calls between-subject shuffling "probably more
  conservative, and suggested for most cases", `R/rMEA_rand.R`.

**bsync departures (stated, not from the source).**
- `keep_roles = TRUE` is the default, because the sign of bsync's lead and
  lag depends on which series is `x`. `FALSE` reproduces the rMEA pair set.
- No `na.omit()`. Removing NA rows shifts the timeline (Invariant 8).
  Missing values pass through to the estimator's own NA policy.
- `generate_surrogate_pseudo()` builds partners for one target dyad. It
  excludes partners shorter than the target's `y`, because the surrogate
  wrappers need `nrow(y_surrogates) == length(y)`. rMEA has no per-dyad form.

## Traces to

- `R/surrogate_generation.R`: the roxygen of `generate_surrogate_pseudo()` and
  `generate_pseudo_dyads()` (M009).
- `vignettes/surrogate-testing.Rmd`: the pseudo-dyad section (M009).
- `tests/testthat/test-surrogate-pseudo.R`: pair counts and crop rule (M009).

## Open questions

- None. The pairing, role, sampling, and crop rules come from the code read
  above — observed 2026-10-05.

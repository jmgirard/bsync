# tschacher2020 — the SUSY segment convention

**Provenance.** Ingested 2026-10-07 by M015 from the SUSY 0.1.0 CRAN source, read through the GitHub CRAN mirror (`gh api repos/cran/SUSY/contents/...`, mirror commit `02590aea0877f09411584827cecaa623e46b6e4f`). Files read: `R/susy.R` (`susy()`, shelf copy `sources/susy-0.1.0-susy.R`), `man/susy.Rd`, and `DESCRIPTION`. The article itself was not read.
Pagination: — (source code. Anchors are file names and the variable names in `susy()`.)
Extraction: verified 2026-10-07 against the SUSY 0.1.0 files named above, read directly — observed 2026-10-07.

**Citation.** Tschacher, W., & Meier, D. (2020). Physiological synchrony in
psychotherapy sessions. *Psychotherapy Research*, 30(5), 558–573.
DOI: 10.1080/10503307.2019.1612114. The DOI is the one that SUSY's
`DESCRIPTION` prints ("'SUSY' works as described in Tschacher & Meier
(2020)"). The volume, issue, and pages match the existing citation in
`vignettes/choosing-parameters.Rmd`. This note did not check them against
the article.

**Role.** This note records the segment convention of SUSY that
`generate_surrogate_segment()` follows, and the places where bsync departs
from it. The roxygen of that function and the surrogate-testing vignette
cite the article. The SUSY mean absolute Z kernel is pinned separately in
`tests/testthat/test-external-oracle.R`.

## Extracted values

- Unit. `segment` and `maxlag` are in seconds. The segment length in rows
  is `range = segment * Hz`, `R/susy.R`.
- Missing values. Each pair is filtered with `complete.cases(a, b)` before
  the segments are cut. `size` is the number of complete rows, `R/susy.R`.
- Segment count. `numberEpochen = round(size/range-0.499999)`, `R/susy.R`.
  If the fraction of size / range is below 0.999999, this equals
  floor(size / range).
- Segment cut. The `eps` loop fills segment `i0` with rows
  `(i0 - 1) * range + 1` to `i0 * range`, counted from row 1, `R/susy.R`.
- Dropped tail. Rows after `numberEpochen * range` are not used. The Rd
  says that the remainder "is not considered", `man/susy.Rd`, Details.
- Segment-size rules. `stop("'segment' must not be smaller than '2 * maxlag'")`
  and `stop("'segment' must not be greater than 'nrow(x)/2'")`, `R/susy.R`.
  The second rule compares `segment` in seconds with `nrow(x)` in rows.
- Surrogate set. With `restrict.surrogates = FALSE` (the default), the
  pseudo pairs are every segment `i` of `a` with every segment `h != i` of
  `b`: k(k - 1) pairs for k segments. The real pairs are `i == h`. With
  `restrict.surrogates = TRUE`, each segment of `a` gets
  `floor(surrogates.total/numberEpochen)` partners, drawn with `runif()`,
  `R/susy.R`.

**bsync departures (stated, not from the source).**
- bsync reorders whole segments of one series into a full-length surrogate,
  with no segment in its original position, so every surrogate wrapper and
  the multiverse can use it. SUSY compares segment pairs and does not build
  a series. The plan gate of M015 chose this form.
- `segment_size` is in samples, from 2 to floor(length(y) / 2). This is the
  SUSY half-length rule stated in rows. The lag rule
  (`segment >= 2 * maxlag`) is not checked, because the generator has no
  lag. In the multiverse the segment size is the window size, and
  `wcc()` already caps `lag_max` at half the window.
- The tail stays in place at the end, so the surrogate keeps `length(y)`
  rows. SUSY drops it.
- `NA` values move with their segments, and the segment count uses
  length(y), not the number of complete cases. If rows are removed, the
  timeline shifts (Invariant 8).

## Traces to

- `R/surrogate_generation.R`: the roxygen and the segment cut of
  `generate_surrogate_segment()` (M015).
- `tests/testthat/test-surrogate-segment.R`: the hand-derived segment
  boundaries from the `eps` loop (M015).
- `vignettes/surrogate-testing.Rmd`: the segment-shuffling section (M015).

## Open questions

- None. The segment count, cut, tail, and size rules come from the code
  read above — observed 2026-10-07.

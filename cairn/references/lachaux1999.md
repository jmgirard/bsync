# lachaux1999 — the phase-locking value (PLV) definition wphase implements

**Provenance.** Ingested 2026-08-29 by M008 from
`cairn/references/sources/lachaux1999.pdf` (gitignored; open-access green copy
retrieved from Europe PMC, `https://europepmc.org/articles/pmc6873296?pdf=render`,
located via the Semantic Scholar API by DOI).
Pagination: journal pages (194–208, printed on the PDF pages).
Extraction: verified 2026-08-29 against the source (pp. 194–198 read directly; the formula and surrogate-test passages quoted below) — observed 2026-08-29.

**Citation.** Lachaux, J.-P., Rodriguez, E., Martinerie, J., & Varela, F. J.
(1999). Measuring phase synchrony in brain signals. *Human Brain Mapping*,
8(4), 194–208. DOI: 10.1002/(SICI)1097-0193(1999)8:4<194::AID-HBM4>3.0.CO;2-C.

**Role.** Primary source for the PLV formula `wphase()` computes and the
surrogate logic `wphase_surrogate()` follows. The oracle tests cite this page;
the roxygen docs name this lineage.

## Extracted values

- PLV definition — "PLV_t = (1/N) |Σ_{n=1}^{N} exp (jθ(t, n))|", "where
  θ(t, n) is the phase difference φ1(t, n) − φ2(t, n)", p. 195. The modulus
  of the mean unit phasor of the phase differences; "If the phase difference
  varies little across the trials, PLV is close to 1; it is close to zero
  otherwise" (p. 196).
- Phase extraction — instantaneous phase at a frequency of interest, obtained
  after band-pass filtering (their step 1) via convolution with a complex
  Gabor wavelet (their step 2), p. 195. Filtering before phase extraction is
  integral to the method: the authors tested skipping the filter ("method A")
  and "relied on method B" (filter first) because method A "occasionally gave
  different results", p. 198.
- Significance — surrogate data ("step 3"), shuffling trials of one electrode;
  "The proportion of surrogate values higher than the original PLV … is
  called the phase-locking statistics (PLS)", 5% criterion, p. 197.
- Known blind spot the authors state — phase differences constant across
  trials with no variability defeat the shuffle null ("except in one
  important case: when the phase values φ1(n) and φ2(n) remains constant
  across the trials"), p. 198.

**bsync adaptation (stated, not from the source).** Lachaux et al. average
the unit phasor across *trials* at a fixed latency. bsync's `wphase()` is a
single-trial continuous-data estimator: it applies the same modulus-of-mean-
phasor formula across *time samples within a sliding window* (at a lag
offset), with phases from `gsignal::hilbert` on the full series — the users'
band-pass responsibility is documented rather than built in (the vignette
candidate row carries narrowband guidance). The surrogate null uses bsync's
matched-null engine (circular shift), not the trial shuffle, which has no
analogue for single-trial data.

## Traces to

- `tests/testthat/test-wphase.R` — closed-form window fixture and analytic
  checks cite the p. 195 formula (file added by M008; lines recorded at
  review).
- `R/wphase.R` roxygen — names the PLV lineage (M008).

## Open questions

- None — the formula and its context are extracted; the adaptation choices
  are recorded above and in M008's plan — observed 2026-08-29.

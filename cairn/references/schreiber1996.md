# schreiber1996 — the IAAFT surrogate algorithm

**Provenance.** Ingested 2026-10-06 by M014 from `cairn/references/sources/schreiber1996.pdf` (gitignored). The file is the arXiv chao-dyn/9909041v1 copy of 30 Sep 1999 (4 pages, 115 KB), fetched from arxiv.org on 2026-10-06 with the user's approval at the M014 plan gate.
Pagination: preprint pages. The arXiv copy has its own page numbers 1 to 4, not the journal pages 635 to 638.
Extraction: verified 2026-10-06 against the source, preprint pp. 1-4 read directly — observed 2026-10-06.

**Citation.** Schreiber, T., & Schmitz, A. (1996). Improved surrogate data for
nonlinearity tests. *Physical Review Letters*, 77(4), 635-638. DOI:
10.1103/PhysRevLett.77.635. The arXiv copy prints the journal reference as
"Phys. Rev. Lett. 77, 635 (1996)". It does not print the issue number, the end
page, or the DOI. Those three follow the APS pattern for this volume and page,
and this note did not read them from the journal.

**Role.** This note is the source for `generate_surrogate_iaaft()`: the
iteration, the stopping rule, the series it returns, and the AR(1) test
process that its calibration test uses.

## Extracted values

- Null hypothesis. A series is consistent with the null "if there exists an
  underlying Gaussian linear stochastic signal {xn} such that sn = h(xn)",
  where h is an "instantaneous, invertible measurement function", p. 1.
- Why AAFT is not enough. AAFT surrogates have the same values as the data
  "by construction, but they do not usually have the same sample power
  spectra", p. 2. AAFT gave "66±5% false rejections" where 5% was expected,
  p. 3.
- Stored targets. "Store a sorted list of the values {sn} and the squared
  amplitudes of the Fourier transform of {sn}", with
  S_k^2 = |sum_{n=0}^{N-1} s_n exp(i 2 pi k n / N)|^2, p. 2.
- Start. "Begin with a random shuffle (without replacement) {s_n^(0)} of the
  data", p. 2. Footnote 10 names a start from an AAFT surrogate as an
  alternative, "Alternatively but not equivalently", p. 4.
- Step 1. Take the Fourier transform of s^(i), replace the squared amplitudes
  by S_k^2, and transform back. "The phases of the complex Fourier components
  are kept", p. 2.
- Step 2. "rank-order the resulting series in order to assume exactly the
  values taken by {sn}", p. 2. The two steps "have to be repeated several
  times", p. 2.
- Stopping. "At each iteration stage we can check the remaining discrepancy
  of the spectrum and iterate until a given accuracy is reached", p. 2. The
  paper also says that "Eventually, the transformation towards the correct
  spectrum will result in a change which is too small to cause a reordering
  of the values. Thus after rescaling, the sequence is not changed", p. 2.
  It advises to "report a failure if this accuracy could not be reached",
  p. 2.
- The returned series. The surrogates "assume the same values (without
  replacement) as the data" and have spectra "practically indistinguishable"
  from the data's, p. 2.
- Discrepancy measure. The relative discrepancy at iteration i is
  sum_k (S_hat_k^(i) - S_hat_k)^2 / sum_k S_hat_k^2, with S_hat_k the
  spectrum smoothed over 21 bins for the figure. "for the generation of
  surrogates no smoothing is performed", pp. 2-3.
- Convergence. The discrepancy falls roughly like 1/i until an N-dependent
  saturation value, p. 3, Fig. 2. The initial discrepancy of a random
  scramble was "0.2±0.01", p. 3.
- Test process. "a first order AR process xn = 0.7xn−1 + ηn, measured through
  sn = x^3_n", with independent Gaussian increments, p. 2.
- Size. With 19 surrogates and a one-sided test, "The correct rejection rate
  for the 95% level of significance is reached after about 7 iterations",
  p. 3, Fig. 4 (a different AR process, sn = xn sqrt|xn|, xn = 0.95xn−1 + ηn).

**bsync choices (stated, not from the source).**
- The stop is the fixed point (the rank-ordering leaves the series
  unchanged), or `max_iter`. The paper's own stop is a spectral accuracy
  target, and it describes the fixed point as where the iteration ends.
- The start is a random shuffle, not an AAFT surrogate.
- The generator returns the rank-ordered series, so each surrogate holds
  exactly the values of `y`.

## Traces to

- `R/surrogate_generation.R`, `generate_surrogate_iaaft()` — the iteration,
  the stop, and the returned series.
- `tests/testthat/test-surrogate-iaaft.R` — the plain reference
  implementation, the discrepancy measure, and the AR(1) cube process of the
  calibration test.

## Open questions

- None.

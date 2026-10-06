# phipson2010 — the add-one (exact) Monte Carlo p-value

**Provenance.** Ingested 2026-10-05 by M011 from `cairn/references/sources/phipson2010.pdf` (gitignored). The file is the arXiv:1603.05766v1 preprint of 18 Mar 2016 (199.6 KB), fetched 2026-10-05.
Pagination: preprint pages. The printed page numbers match the PDF page numbers.
Extraction: verified 2026-10-05 against the source, PDF pp. 1, 3–7 read directly — observed 2026-10-05.

**Citation.** Phipson, B., & Smyth, G. K. (2010). Permutation p-values should
never be zero: calculating exact p-values when permutations are randomly drawn.
*Statistical Applications in Genetics and Molecular Biology*, 9(1), Article 39.
DOI: 10.2202/1544-6115.1585. This is the form the preprint's p. 1 box asks
readers to cite. The preprint title capitalizes it as "Permutation P-values
Should Never Be Zero: Calculating Exact P-values When Permutations Are Randomly
Drawn", and its p. 1 says "Published 31 October 2010; corrected 9 February 2011".

**Role.** This note is the source for the surrogate p-value formula that every
`*_surrogate()` wrapper returns (D-001). It also gives the reason for the
p <= .05 significance rule (D-002): test size is defined as P(p̂ ≤ α) on p. 4,
and the add-one form gives "a test with the correct size" on p. 6. The
few-surrogates warning in `synchrony_multiverse()` fires below 19 surrogates,
because the smallest add-one p-value is 1 / (n + 1) (corrected M012: was
`n_surrogates < 20`).

## Extracted values

- The b / m estimator. "suppose that p̂ = B/m where B is the number of Monte
  Carlo test statistics greater than or equal to t", p. 3 (§2).
- Size of b / m. Under the null, p̂ is discrete uniform on 0, 1/m, …, 1, with
  P(p̂ = b/m) = 1 / (m + 1), and P(p̂ ≤ α) = (⌊mα⌋ + 1) / (m + 1), p. 4 (§2).
  "the true type I error rate P(p̂ ≤ α) exceeds α for most practical values of
  m and α", p. 4. P(p̂ ≤ α) "is never less than 1/(m+1)", p. 4.
- The exact Monte Carlo p-value. B counts the m simulated statistics with
  t_sim ≥ t_obs, and under the null B is discrete uniform on 0, …, m. "Hence
  the exact Monte Carlo p-value is p_u = P(B ≤ b) = (b+1)/(m+1)", p. 6 (§4).
  Its positive bias gives "a test with the correct size", p. 6.
- Exhaustive enumeration. If all m_t possible distinct permutations are
  evaluated, "The exhaustive permutation p-value is p_t = (b_t+1)/(m_t+1)",
  p. 7 (§5). The exhaustive pseudo-dyad null, where every other dyad is a
  partner, is this case.
- Sampling with replacement. If permutations are drawn with replacement, the
  exact p-value "is now slightly less than (b+1)/(m+1)", p. 7 (§6.1). For that
  case the add-one form is conservative.

## Traces to

- `R/surrogate_analysis.R`, the p-value line in `wcc_surrogate()`,
  `wdtw_surrogate()`, `wgranger_surrogate()` (both directions), and
  `wphase_surrogate()`. These compute (b + 1) / (n + 1), p. 6.
- `R/surrogate_analysis.R`, `is_significant()` (rejects at p <= .05, pp. 4, 6)
  and `significance_reachable()` (the smallest add-one p-value, 1 / (n + 1)).
- `R/multiverse.R`, `warn_few_surrogates()`. It uses the smallest add-one
  p-value, 1 / (n + 1).
- `tests/testthat/test-surrogate.R`, the tests "all four wrappers return
  (b + 1) / (n + 1) when b = 0" and "... when 0 < b < n" (closed-form oracle,
  p. 6), and the seeded pin "AC4: seeded surrogate p-values are reproducible".
- `tests/testthat/test-surrogate-pseudo.R`, the test "the add-one pseudo-dyad
  p-value keeps its size under the null" (simulation-coverage oracle,
  pp. 4, 6).
- `cairn/DESIGN.md` §6, step 3 and the Invariant 2 paragraph.
- `cairn/DECISIONS.md` D-001.

## Open questions

- None — observed 2026-10-05.

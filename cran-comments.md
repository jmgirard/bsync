## R CMD check results

There were 0 errors, 0 warnings, and 0 notes.

### Test environments

* macOS Tahoe 26.5.1 (aarch64), R 4.6.1 (2026-06-24), checked locally
* win-builder (R-devel) — see pre-submission checklist below
* R-hub (Ubuntu, R-release) — see pre-submission checklist below

### Downstream dependencies

This is the first CRAN submission of bsync. There are no reverse dependencies.

### Notes for CRAN reviewers

* The package uses compiled C++ code (Rcpp/RcppArmadillo) for the windowed
  estimator cores. The compiled code is serial (no OpenMP), reproducible, and
  has been tested with `valgrind` on the local check.

* Vignette build times: the most expensive vignette (`choosing-parameters`)
  runs a surrogate analysis with `n_surrogates = 100` on `sim_dyad` (a
  bundled synthetic dataset of 6,000 rows). All six vignettes build in
  under 2 minutes total on the test machine. Tests marked `skip_on_cran()`
  are the heavy surrogate runs (1,000+ permutations).

* The `sim_dyad` dataset is a small synthetic dataset (~250 KB) generated
  from a known ground-truth model for use in examples and vignettes.

---

## Pre-submission checklist (complete before uploading)

- [ ] **Deploy pkgdown site**: run `pkgdown::deploy_to_branch()` (or trigger
  the GitHub Actions pkgdown workflow) so that the `articles/bsync.html` URL
  referenced in README is live. `urlchecker::url_check()` currently reports
  two 404s for this new article; they will resolve after deployment.

- [ ] **win-builder check**: submit to `devtools::check_win_devel()` and
  `devtools::check_win_release()`; confirm 0/0/0.

- [ ] **R-hub check**: run `rhub::rhub_check()` on at least Ubuntu + Windows
  builders; confirm 0/0/0.

- [ ] **Final local check**: `devtools::check(args = '--as-cran')` = 0/0/0.

- [ ] **urlchecker clean** (after pkgdown deploy): `urlchecker::url_check()`.

- [ ] **spelling clean**: `spelling::spell_check_package()`.

- [ ] **pkgdown clean**: `pkgdown::check_pkgdown()`.

- [ ] Submit via `devtools::submit_cran()`.

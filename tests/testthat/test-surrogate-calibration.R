library(testthat)

# Size (Type-I) calibration of surrogate generators (DESIGN.md section 13,
# layer 3). Under an independent pair, a valid null rejects at p <= .05 about
# 5% of the time.
#
# Oracle records:
#   iaaft-calibration (simulation-coverage): 400 independent pairs of length
#     512. Each partner is an AR(1) series x_n = 0.7 x_{n-1} + eta_n observed
#     as s_n = x_n^3 (Schreiber & Schmitz, 1996, p. 2; source note
#     cairn/references/schreiber1996.md). wcc_surrogate() with window 32,
#     lag 4, increment 16, and 19 IAAFT surrogates per pair. Asserted below.

ar_cube_pair <- function(n, phi = 0.7) {
  draw <- function() {
    x <- as.numeric(stats::filter(stats::rnorm(n), phi, method = "recursive"))
    x^3
  }
  list(x = draw(), y = draw())
}

iaaft_null_p <- function(match, n_pairs = 400, n = 512, seed = 20261006) {
  set.seed(seed)
  vapply(
    seq_len(n_pairs),
    function(i) {
      d <- ar_cube_pair(n)
      surr <- generate_surrogate_iaaft(d$y, n_surrogates = 19, match = match)
      wcc_surrogate(
        d$x,
        d$y,
        surr,
        window_size = 32,
        lag_max = 4,
        window_increment = 16
      )$p_value
    },
    numeric(1)
  )
}

test_that("IAAFT spectrum form keeps WCC test size; values form does not", {
  skip_on_cran()
  p_spectrum <- iaaft_null_p("spectrum")
  expect_length(p_spectrum, 400)
  expect_lte(mean(p_spectrum <= 0.05), 0.08)
  # A null that collapsed toward y (p near 1) would also pass the bound.
  expect_gte(mean(p_spectrum <= 0.05), 0.01)

  # The generator draws the same random numbers for both forms, so the first
  # 200 pairs here are the first 200 pairs above.
  p_values <- iaaft_null_p("values", n_pairs = 200)
  rate_values <- mean(p_values <= 0.05)
  expect_gt(rate_values, 0.05)
  expect_gt(rate_values, mean(p_spectrum[1:200] <= 0.05))
})

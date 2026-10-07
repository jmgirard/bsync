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
#   segment-calibration (simulation-coverage): 1000 independent pairs of
#     AR(1) series x_n = 0.7 x_{n-1} + eta_n per run, at 5 segments (length
#     100, the enumeration path) and 12 segments (length 240, the rejection
#     path). wcc_surrogate() with window 16, lag 2, and increment 20, and 19
#     segment surrogates of segment_size 20, so that every lagged window lies
#     inside one segment. Design and bounds: RR01 Q6
#     (cairn/reviews/archive/RR01-segment-surrogate-size.md); order rule:
#     D-004. A control with orders restricted to no segment in place shows
#     the test can fail. Asserted below.

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

ar1_pair <- function(n, phi = 0.7) {
  draw <- function() {
    as.numeric(stats::filter(stats::rnorm(n), phi, method = "recursive"))
  }
  list(x = draw(), y = draw())
}

# Surrogate matrix of y from a matrix of segment orders, one order per row.
segments_from_orders <- function(y, orders, segment_size) {
  k <- ncol(orders)
  tail_rows <- seq_len(length(y) - k * segment_size) + k * segment_size
  apply(orders, 1, function(o) {
    y[c(
      rep((o - 1) * segment_size, each = segment_size) +
        seq_len(segment_size),
      tail_rows
    )]
  })
}

segment_null_p <- function(n, seed, n_pairs = 1000, no_fixed = FALSE) {
  set.seed(seed)
  vapply(
    seq_len(n_pairs),
    function(i) {
      d <- ar1_pair(n)
      surr <- if (no_fixed) {
        # Planted defect: draw 19 of the orders with no segment in place.
        orders <- .segment_orders(n %/% 20)
        orders <- orders[apply(orders, 1, function(o) all(o != seq_along(o))), ]
        segments_from_orders(d$y, orders[sample.int(nrow(orders), 19), ], 20)
      } else {
        generate_surrogate_segment(d$y, segment_size = 20, n_surrogates = 19)
      }
      wcc_surrogate(
        d$x,
        d$y,
        surr,
        window_size = 16,
        lag_max = 2,
        window_increment = 20,
        statistic = "mean_abs_z"
      )$p_value
    },
    numeric(1)
  )
}

test_that("segment surrogates keep WCC test size with aligned windows", {
  skip_on_cran()
  p_5 <- segment_null_p(100, seed = 20261007)
  expect_length(p_5, 1000)
  expect_gte(mean(p_5 <= 0.05), 0.03)
  expect_lte(mean(p_5 <= 0.05), 0.07)

  p_12 <- segment_null_p(240, seed = 20261008)
  expect_length(p_12, 1000)
  expect_gte(mean(p_12 <= 0.05), 0.03)
  expect_lte(mean(p_12 <= 0.05), 0.07)

  # Control: 44 of the 119 orders of 5 segments keep no segment in place.
  # Drawing only from them rejects a true null too often (RR01 Q1).
  p_control <- segment_null_p(100, seed = 20261007, no_fixed = TRUE)
  expect_gt(mean(p_control <= 0.05), 0.07)
})

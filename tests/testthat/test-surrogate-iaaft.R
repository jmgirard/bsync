library(testthat)

# IAAFT surrogates. Source: Schreiber & Schmitz (1996), Phys. Rev. Lett. 77,
# 635; source note cairn/references/schreiber1996.md.
#
# Oracle records (DESIGN.md section 13):
#   iaaft-closed-form (closed-form): iaaft_reference() below, a plain
#     reimplementation of the p. 2 iteration. It uses the explicit DFT sum of
#     the paper's S_k formula (sign +i, so it also checks that the result does
#     not depend on the transform's sign convention), explicit loops, the
#     fixed-point stop, and the generator's random-number order: one
#     sample.int(n) shuffle per column, column 1 first. Asserted in the
#     "matches the plain reference" test.
#   iaaft-invariants (invariant): exact value distribution, fixed point after
#     convergence, smaller spectral discrepancy than a shuffle (pp. 2-3).
#   iaaft-calibration (simulation-coverage): test-surrogate-calibration.R.

# --- Plain reference (schreiber1996, p. 2) ---------------------------------

# Explicit DFT with the paper's sign: X_k = sum_n s_n exp(+i 2 pi k n / N).
dft_plain <- function(s) {
  n <- length(s)
  out <- complex(n)
  for (k in 0:(n - 1)) {
    acc <- 0 + 0i
    for (m in 0:(n - 1)) {
      acc <- acc + s[m + 1] * exp(1i * 2 * pi * k * m / n)
    }
    out[k + 1] <- acc
  }
  out
}

# Inverse of dft_plain: s_n = (1 / N) sum_k X_k exp(-i 2 pi k n / N).
idft_plain <- function(x) {
  n <- length(x)
  out <- complex(n)
  for (m in 0:(n - 1)) {
    acc <- 0 + 0i
    for (k in 0:(n - 1)) {
      acc <- acc + x[k + 1] * exp(-1i * 2 * pi * k * m / n)
    }
    out[m + 1] <- acc / n
  }
  Re(out)
}

# One iteration: impose the data's Fourier amplitudes, keep the phases, then
# rank-order so that the series takes exactly the data's values.
iaaft_step <- function(s, y) {
  target_amp <- Mod(dft_plain(y))
  spec <- dft_plain(s)
  adjusted <- idft_plain(target_amp * exp(1i * Arg(spec)))
  y_sorted <- sort(y)
  out <- numeric(length(s))
  ranks <- rank(adjusted, ties.method = "first")
  for (m in seq_along(s)) {
    out[m] <- y_sorted[ranks[m]]
  }
  out
}

iaaft_reference <- function(y, n_surrogates, max_iter) {
  n <- length(y)
  out <- matrix(0, nrow = n, ncol = n_surrogates)
  starts <- vector("list", n_surrogates)
  for (j in seq_len(n_surrogates)) {
    starts[[j]] <- y[sample.int(n)]
  }
  for (j in seq_len(n_surrogates)) {
    s <- starts[[j]]
    for (i in seq_len(max_iter)) {
      s_new <- iaaft_step(s, y)
      done <- all(s_new == s)
      s <- s_new
      if (done) {
        break
      }
    }
    out[, j] <- s
  }
  out
}

# Unsmoothed relative spectral discrepancy (schreiber1996, pp. 2-3).
spectral_discrepancy <- function(s, y) {
  s2 <- Mod(dft_plain(s))^2
  y2 <- Mod(dft_plain(y))^2
  sum((s2 - y2)^2) / sum(y2^2)
}

# A correlated, skewed test series: AR(1) observed through a cube.
ar_cube <- function(n, phi = 0.7) {
  x <- as.numeric(stats::filter(stats::rnorm(n), phi, method = "recursive"))
  x^3
}

# --- AC1: shape and exact value distribution -------------------------------

test_that("IAAFT surrogates have the data's exact values (even, odd, ties)", {
  set.seed(101)
  inputs <- list(
    even = ar_cube(64),
    odd = ar_cube(63),
    ties = round(ar_cube(50), 1)
  )
  expect_gt(anyDuplicated(inputs$ties), 0)
  for (y in inputs) {
    surr <- generate_surrogate_iaaft(y, n_surrogates = 7)
    expect_true(is.matrix(surr))
    expect_true(is.numeric(surr))
    expect_identical(dim(surr), c(length(y), 7L))
    for (j in seq_len(ncol(surr))) {
      expect_identical(sort(surr[, j]), sort(y))
    }
  }
})

test_that("IAAFT accepts integer input and returns doubles", {
  set.seed(102)
  y <- sample.int(1000L, 40)
  surr <- generate_surrogate_iaaft(y, n_surrogates = 3)
  expect_type(surr, "double")
  expect_identical(sort(surr[, 1]), sort(as.double(y)))
})

# --- AC2: closed-form oracle -----------------------------------------------

test_that("IAAFT matches the plain reference (schreiber1996, p. 2)", {
  for (case in list(list(seed = 201, n = 32), list(seed = 202, n = 25))) {
    set.seed(case$seed)
    y <- ar_cube(case$n)

    set.seed(case$seed + 1000)
    got <- generate_surrogate_iaaft(y, n_surrogates = 4)
    set.seed(case$seed + 1000)
    want <- iaaft_reference(y, n_surrogates = 4, max_iter = 1000)

    expect_equal(got, want)
  }
})

# --- AC3: fixed point and spectrum -----------------------------------------

test_that("converged IAAFT columns are fixed points with a closer spectrum", {
  set.seed(301)
  y <- ar_cube(32)
  surr <- expect_no_warning(
    generate_surrogate_iaaft(y, n_surrogates = 5, max_iter = 1000)
  )
  shuffle <- y[sample.int(length(y))]
  d_shuffle <- spectral_discrepancy(shuffle, y)
  for (j in seq_len(ncol(surr))) {
    expect_identical(iaaft_step(surr[, j], y), surr[, j])
    expect_lt(spectral_discrepancy(surr[, j], y), d_shuffle)
  }
})

# --- AC5: warning and errors -----------------------------------------------

test_that("IAAFT warns about unconverged columns and still returns them", {
  set.seed(501)
  y <- ar_cube(200)
  expect_warning(
    surr <- generate_surrogate_iaaft(y, n_surrogates = 6, max_iter = 1),
    "6 of 6 surrogates did not converge"
  )
  expect_identical(dim(surr), c(200L, 6L))
  expect_identical(sort(surr[, 6]), sort(y))
})

test_that("IAAFT rejects invalid input", {
  y <- c(1, 4, 2, 8, 5, 7)
  expect_error(
    generate_surrogate_iaaft(c(y, NA)),
    "must not contain missing values"
  )
  expect_error(
    generate_surrogate_iaaft(c(y, Inf)),
    "must contain only finite values"
  )
  expect_error(
    generate_surrogate_iaaft(as.character(y)),
    "must be a numeric vector"
  )
  expect_error(generate_surrogate_iaaft(c(1, 2)), "at least 3 values")
  expect_identical(dim(generate_surrogate_iaaft(c(1, 3, 2), 2)), c(3L, 2L))
  expect_error(
    generate_surrogate_iaaft(y, n_surrogates = 0),
    "n_surrogates.*single positive integer"
  )
  expect_error(
    generate_surrogate_iaaft(y, n_surrogates = 2.5),
    "n_surrogates.*single positive integer"
  )
  expect_error(
    generate_surrogate_iaaft(y, max_iter = 0),
    "max_iter.*single positive integer"
  )
  expect_error(
    generate_surrogate_iaaft(y, max_iter = c(10, 20)),
    "max_iter.*single positive integer"
  )
})

test_that("the IAAFT missing-value error points to impute_ts_gaps()", {
  expect_error(
    generate_surrogate_iaaft(c(1, 4, NA, 2, 8)),
    "impute_ts_gaps"
  )
})

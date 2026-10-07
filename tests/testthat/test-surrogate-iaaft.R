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
#     sample.int(n) shuffle per column, column 1 first. It returns both
#     forms: the rank-ordered series and the last spectrum-adjusted series.
#     Asserted in the "matches the plain reference" test.
#   iaaft-invariants (invariant): exact values (match = "values") or exact
#     Fourier amplitudes (match = "spectrum"), fixed point after convergence,
#     smaller spectral discrepancy than a shuffle (pp. 2-3), both forms from
#     one iteration, spectrum-form values closer to y than phase surrogates.
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

# Step 1: impose the data's Fourier amplitudes and keep the phases.
iaaft_adjust <- function(s, y) {
  target_amp <- Mod(dft_plain(y))
  spec <- dft_plain(s)
  idft_plain(target_amp * exp(1i * Arg(spec)))
}

# Step 2: rank-order so that the series takes exactly the data's values.
iaaft_rank <- function(adjusted, y) {
  y_sorted <- sort(y)
  out <- numeric(length(adjusted))
  ranks <- rank(adjusted, ties.method = "first")
  for (m in seq_along(adjusted)) {
    out[m] <- y_sorted[ranks[m]]
  }
  out
}

# One full iteration.
iaaft_step <- function(s, y) {
  iaaft_rank(iaaft_adjust(s, y), y)
}

# Returns both forms: `values` (the rank-ordered series) and `spectrum` (the
# last spectrum-adjusted series, before the rank step).
iaaft_reference <- function(y, n_surrogates, max_iter) {
  n <- length(y)
  values <- matrix(0, nrow = n, ncol = n_surrogates)
  spectrum <- matrix(0, nrow = n, ncol = n_surrogates)
  starts <- vector("list", n_surrogates)
  for (j in seq_len(n_surrogates)) {
    starts[[j]] <- y[sample.int(n)]
  }
  for (j in seq_len(n_surrogates)) {
    s <- starts[[j]]
    for (i in seq_len(max_iter)) {
      adjusted <- iaaft_adjust(s, y)
      s_new <- iaaft_rank(adjusted, y)
      done <- all(s_new == s)
      s <- s_new
      if (done) {
        break
      }
    }
    values[, j] <- s
    spectrum[, j] <- adjusted
  }
  list(values = values, spectrum = spectrum)
}

# Relative spectral discrepancy (schreiber1996, p. 3) without the paper's
# 21-bin smoothing: sum_k (S_k(s) - S_k(y))^2 / sum_k S_k(y)^2, where S_k is
# the square root of the power at bin k (the Fourier amplitude).
spectral_discrepancy <- function(s, y) {
  s_amp <- Mod(dft_plain(s))
  y_amp <- Mod(dft_plain(y))
  sum((s_amp - y_amp)^2) / sum(y_amp^2)
}

# A correlated, skewed test series: AR(1) observed through a cube.
ar_cube <- function(n, phi = 0.7) {
  x <- as.numeric(stats::filter(stats::rnorm(n), phi, method = "recursive"))
  x^3
}

# --- AC1: shape and exact value distribution -------------------------------

iaaft_inputs <- local({
  set.seed(101)
  list(
    even = ar_cube(64),
    odd = ar_cube(63),
    ties = round(ar_cube(50), 1)
  )
})

test_that("IAAFT with match = 'values' keeps the exact values", {
  expect_gt(anyDuplicated(iaaft_inputs$ties), 0)
  for (y in iaaft_inputs) {
    surr <- generate_surrogate_iaaft(y, n_surrogates = 7, match = "values")
    expect_true(is.matrix(surr))
    expect_true(is.numeric(surr))
    expect_identical(dim(surr), c(length(y), 7L))
    for (j in seq_len(ncol(surr))) {
      expect_identical(sort(surr[, j]), sort(y))
    }
  }
})

test_that("IAAFT with match = 'spectrum' (default) keeps the exact spectrum", {
  for (y in iaaft_inputs) {
    surr <- generate_surrogate_iaaft(y, n_surrogates = 7)
    expect_true(is.matrix(surr))
    expect_true(is.numeric(surr))
    expect_identical(dim(surr), c(length(y), 7L))
    for (j in seq_len(ncol(surr))) {
      expect_equal(Mod(stats::fft(surr[, j])), Mod(stats::fft(y)))
    }
  }
})

test_that("IAAFT accepts integer input and returns doubles", {
  set.seed(102)
  y <- sample.int(1000L, 40)
  surr <- generate_surrogate_iaaft(y, n_surrogates = 3, match = "values")
  expect_type(surr, "double")
  expect_identical(sort(surr[, 1]), sort(as.double(y)))
})

# --- AC2: closed-form oracle -----------------------------------------------

test_that("IAAFT matches the plain reference (schreiber1996, p. 2)", {
  for (case in list(list(seed = 201, n = 32), list(seed = 202, n = 25))) {
    set.seed(case$seed)
    y <- ar_cube(case$n)

    set.seed(case$seed + 1000)
    want <- iaaft_reference(y, n_surrogates = 4, max_iter = 1000)
    for (form in c("values", "spectrum")) {
      set.seed(case$seed + 1000)
      got <- generate_surrogate_iaaft(y, n_surrogates = 4, match = form)
      expect_equal(got, want[[form]])
    }
  }
})

# --- AC3: fixed point, spectrum, and closeness of values -------------------

test_that("converged IAAFT columns are fixed points with a closer spectrum", {
  set.seed(301)
  y <- ar_cube(32)
  set.seed(302)
  surr <- expect_no_warning(
    generate_surrogate_iaaft(y, n_surrogates = 5, match = "values")
  )
  set.seed(302)
  surr_spec <- generate_surrogate_iaaft(y, n_surrogates = 5)
  shuffle <- y[sample.int(length(y))]
  d_shuffle <- spectral_discrepancy(shuffle, y)
  for (j in seq_len(ncol(surr))) {
    expect_identical(iaaft_step(surr[, j], y), surr[, j])
    expect_lt(spectral_discrepancy(surr[, j], y), d_shuffle)
    # Both forms come from the same final iteration.
    expect_identical(
      sort(y)[rank(surr_spec[, j], ties.method = "first")],
      surr[, j]
    )
  }
})

test_that("spectrum-form values are closer to y than phase surrogates are", {
  set.seed(303)
  y <- ar_cube(128)
  surr <- generate_surrogate_iaaft(y, n_surrogates = 5)
  phase <- generate_surrogate_phase(y, n_surrogates = 20)
  value_dist <- function(s) mean(abs(sort(s) - sort(y)))
  d_phase <- mean(apply(phase, 2, value_dist))
  for (j in seq_len(ncol(surr))) {
    expect_lt(value_dist(surr[, j]), d_phase)
  }
})

# --- AC5: warning and errors -----------------------------------------------

test_that("IAAFT warns about unconverged columns and still returns them", {
  set.seed(501)
  y <- ar_cube(200)
  hints <- c(
    spectrum = "for values closer to the values of",
    values = "for a spectrum closer to the spectrum of"
  )
  for (form in names(hints)) {
    w <- expect_warning(
      surr <- generate_surrogate_iaaft(
        y,
        n_surrogates = 6,
        max_iter = 1,
        match = form
      ),
      "6 of 6 surrogates did not converge"
    )
    expect_match(conditionMessage(w), hints[[form]])
    expect_identical(dim(surr), c(200L, 6L))
  }
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
  expect_error(
    generate_surrogate_iaaft(y, match = "phase"),
    "match.*must be one of"
  )
  expect_error(
    generate_surrogate_iaaft(y, match = c("values", "spectrum")),
    "match.*must be a single value"
  )
})

test_that("the IAAFT warning is classed and pluralized", {
  set.seed(502)
  y <- ar_cube(200)
  w <- expect_warning(
    generate_surrogate_iaaft(y, n_surrogates = 1, max_iter = 1),
    class = "bsync_iaaft_unconverged"
  )
  expect_match(conditionMessage(w), "1 of 1 surrogate did not converge")
  expect_identical(w$n_unconverged, 1L)
  expect_identical(w$max_iter, 1)
})

test_that("IAAFT rejects a series whose Fourier transform overflows", {
  expect_error(
    generate_surrogate_iaaft(c(1e308, -1e308, 1e308, 2, 3)),
    "Fourier transform of .*y.* overflows"
  )
})

test_that("the IAAFT missing-value error points to impute_ts_gaps()", {
  expect_error(
    generate_surrogate_iaaft(c(1, 4, NA, 2, 8)),
    "impute_ts_gaps"
  )
})

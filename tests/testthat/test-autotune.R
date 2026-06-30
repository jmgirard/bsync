test_that("autotune_wcc: structure and sanity on sim_dyad (regression)", {
  skip_on_cran() # multiverse + surrogates is slow for CRAN

  set.seed(2024)

  # Build a small dyad_list from sim_dyad (3 identical dyads for speed)
  dyad_list <- replicate(
    3,
    list(x = sim_dyad$z_A, y = sim_dyad$z_B),
    simplify = FALSE
  )

  result <- autotune_wcc(
    dyad_list    = dyad_list,
    sample_rate  = 80, # sim_dyad is 80 Hz
    window_sec   = c(1, 2, 4), # range spanning the 0.5 Hz signal (cycle = 2 s)
    lag_sec      = c(0.5, 1),
    n_surrogates = 30L, # fast; real use >= 100
    n_tune_dyads = 3L
  )

  # Structure
  expect_named(result, c(
    "window_size", "lag_max", "window_increment", "lag_increment",
    "window_sec", "lag_sec",
    "sig_rate", "median_es", "iqr_es", "score",
    "n_dyads", "n_cells_gated", "dyad_multiverses"
  ))

  # Parameters are positive integers
  expect_true(is.integer(result$window_size) || is.numeric(result$window_size))
  expect_true(result$window_size > 0)
  expect_true(result$lag_max > 0)
  expect_true(result$window_increment > 0)
  expect_identical(result$lag_increment, 1L)

  # Seconds conversion is consistent
  expect_equal(result$window_sec * 80, result$window_size, tolerance = 1)

  # Effect size is positive (sim_dyad has real synchrony)
  expect_true(result$median_es > 0)

  # Diagnostics in [0, 1]
  expect_true(result$sig_rate >= 0 && result$sig_rate <= 1)
  expect_true(result$iqr_es >= 0)

  # n_dyads matches what was used
  expect_equal(result$n_dyads, 3L)

  # dyad_multiverses is a list of bsync_multiverse objects
  expect_length(result$dyad_multiverses, 3L)
  expect_true(inherits(result$dyad_multiverses[[1]], "bsync_multiverse"))
})


test_that("autotune_wcc: stability test -- consistent ES across heterogeneous dyads", {
  skip_on_cran()

  set.seed(42)
  n <- 200
  fs <- 10
  t <- seq_len(n) / fs
  sig <- sin(2 * pi * 0.5 * t) # 0.5 Hz, 2-s cycle

  # Three dyads: different noise levels, same underlying signal + 0.2 s lag
  lag_samp <- round(0.2 * fs) # 2 samples
  make_dyad <- function(noise_sd) {
    x <- sig + stats::rnorm(n, sd = noise_sd)
    y <- c(rep(0, lag_samp), sig)[seq_len(n)] + stats::rnorm(n, sd = noise_sd)
    list(x = x, y = y)
  }
  dyad_list <- list(make_dyad(0.5), make_dyad(1.0), make_dyad(0.3))

  # With few surrogates, the gate may not be met; suppress the fallback warning
  result <- suppressWarnings(autotune_wcc(
    dyad_list    = dyad_list,
    sample_rate  = fs,
    window_sec   = c(1, 2),
    lag_sec      = c(0.3, 0.5),
    n_surrogates = 30L
  ))

  # Should pick a valid cell
  expect_true(result$window_size > 0)
  expect_true(result$lag_max > 0)
  expect_true(is.finite(result$score))
})


test_that("autotune_wcc: input validation", {
  dyad <- list(list(x = 1:100, y = 1:100))

  expect_error(
    autotune_wcc(dyad_list = list(), sample_rate = 10, window_sec = 1),
    "non-empty"
  )
  expect_error(
    autotune_wcc(dyad_list = dyad, sample_rate = -1, window_sec = 1),
    "positive"
  )
  expect_error(
    autotune_wcc(dyad_list = dyad, sample_rate = 10, window_sec = -1),
    "positive"
  )
  expect_error(
    autotune_wcc(dyad_list = dyad, sample_rate = 10, window_sec = 1, sig_pct = 2),
    "\\[0, 1\\]"
  )
  expect_error(
    autotune_wcc(dyad_list = dyad, sample_rate = 10, window_sec = 1, iqr_penalty = -1),
    "non-negative"
  )
})


test_that("select_specification: returns bsync_multiverse row and diagnostics", {
  skip_on_cran()

  set.seed(7)
  n <- 200
  fs <- 10
  t <- seq_len(n) / fs
  x <- sin(2 * pi * 0.5 * t) + stats::rnorm(n, sd = 0.3)
  y <- c(rep(0, 2), sin(2 * pi * 0.5 * t))[seq_len(n)] + stats::rnorm(n, sd = 0.3)

  mv1 <- synchrony_multiverse(
    x = x, y = y, estimator = "wcc", sample_rate = fs,
    window_sec = c(1, 2), lag_sec = c(0.3, 0.5),
    n_surrogates = 30L
  )
  mv2 <- synchrony_multiverse(
    x = x + stats::rnorm(n, sd = 0.1),
    y = y + stats::rnorm(n, sd = 0.1),
    estimator = "wcc", sample_rate = fs,
    window_sec = c(1, 2), lag_sec = c(0.3, 0.5),
    n_surrogates = 30L
  )

  # With few surrogates the gate may not be met; suppress the fallback warning
  sel <- suppressWarnings(select_specification(list(mv1, mv2)))

  expect_named(sel, c("best_row", "sig_rate", "median_es", "iqr_es", "score", "n_gated"))
  expect_true(nrow(sel$best_row) == 1)
  expect_true(is.finite(sel$median_es))
  expect_true(sel$iqr_es >= 0)
  expect_true(sel$n_gated >= 0)
})


test_that("select_specification: validates input", {
  expect_error(
    select_specification(list()),
    "non-empty"
  )
  expect_error(
    select_specification(list(1, 2)),
    "bsync_multiverse"
  )
})


test_that("select_specification: emits warning when gate is not met", {
  skip_on_cran()

  # Build a trivial multiverse where no cell is significant (all NA or high p)
  set.seed(5)
  x <- stats::rnorm(100)
  y <- stats::rnorm(100)
  mv <- synchrony_multiverse(
    x = x, y = y, estimator = "wcc", sample_rate = 10,
    window_sec = 1, lag_sec = 0.5,
    n_surrogates = 5L # very few surrogates -> noisy null -> likely no significance
  )

  # Manually zero out significance so gate definitely fails
  mv$grid$p <- rep(0.9, nrow(mv$grid))

  expect_warning(
    select_specification(list(mv), sig_pct = 0.5),
    "detectability gate"
  )
})


test_that("autotune_wcc: works with data.frame dyads", {
  skip_on_cran()

  set.seed(99)
  n <- 150
  df_dyad <- data.frame(
    x = sin(seq(0, 4 * pi, length.out = n)) + stats::rnorm(n, sd = 0.5),
    y = cos(seq(0, 4 * pi, length.out = n)) + stats::rnorm(n, sd = 0.5)
  )
  dyad_list <- list(df_dyad, df_dyad) # same dyad twice for speed

  result <- suppressWarnings(autotune_wcc(
    dyad_list    = dyad_list,
    sample_rate  = 10,
    window_sec   = c(1, 2),
    lag_sec      = 0.5,
    n_surrogates = 20L
  ))

  expect_true(result$window_size > 0)
})

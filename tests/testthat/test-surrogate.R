library(testthat)

# =========================================================================
# --- 1. GENERATOR TESTS --------------------------------------------------
# =========================================================================

test_that("Surrogate generators return matrices of correct dimensions", {
  y <- 1:50
  n_surr <- 15

  surr_circ <- generate_surrogate_circular(y, n_surrogates = n_surr)
  expect_true(is.matrix(surr_circ))
  expect_equal(dim(surr_circ), c(50, 15))

  surr_phase <- generate_surrogate_phase(y, n_surrogates = n_surr)
  expect_true(is.matrix(surr_phase))
  expect_equal(dim(surr_phase), c(50, 15))
})

test_that("Circular shift strictly preserves data distribution", {
  y <- rnorm(50)
  surr_circ <- generate_surrogate_circular(y, n_surrogates = 5)

  expect_equal(sum(surr_circ[, 1]), sum(y))
  expect_equal(mean(surr_circ[, 3]), mean(y))
  expect_true(surr_circ[1, 1] %in% y)
})

test_that("Phase randomization preserves mean and variance", {
  y <- rnorm(100, mean = 5, sd = 2)
  surr_phase <- generate_surrogate_phase(y, n_surrogates = 5)

  expect_equal(mean(surr_phase[, 1]), mean(y), tolerance = 1e-10)
  expect_equal(var(surr_phase[, 1]), var(y), tolerance = 1e-10)
})

test_that("Phase randomization preserves the sign of a negative mean", {
  # Regression: taking Mod() of the DC term used to flip the sign of the mean
  # for negative-mean signals. DC (and Nyquist) must be preserved exactly.
  y <- rnorm(100, mean = -5, sd = 2)
  surr_phase <- generate_surrogate_phase(y, n_surrogates = 5)

  expect_equal(colMeans(surr_phase), rep(mean(y), 5), tolerance = 1e-10)
  expect_equal(var(surr_phase[, 1]), var(y), tolerance = 1e-10)
})

test_that("Phase randomization handles odd-length series natively", {
  # Regression: odd-length input used to abort unless trim_odd = TRUE, which
  # broke the surrogate-matrix length contract in synchrony_multiverse().
  y <- rnorm(101, mean = -2, sd = 1.5)
  surr_phase <- generate_surrogate_phase(y, n_surrogates = 7)

  expect_equal(dim(surr_phase), c(101L, 7L))
  # Output is real and preserves the amplitude spectrum + mean/variance.
  expect_equal(max(abs(Im(surr_phase))), 0)
  expect_equal(Mod(fft(surr_phase[, 1])), Mod(fft(y)), tolerance = 1e-8)
  expect_equal(mean(surr_phase[, 1]), mean(y), tolerance = 1e-10)
  expect_equal(var(surr_phase[, 1]), var(y), tolerance = 1e-10)

  # Legacy trim_odd path still drops the final observation.
  expect_message(
    trimmed <- generate_surrogate_phase(y, n_surrogates = 3, trim_odd = TRUE),
    "Trimming the final observation"
  )
  expect_equal(nrow(trimmed), 100L)
})

test_that("synchrony_multiverse runs on odd-length series with phase surrogates", {
  # Regression: the default surrogate_method = "phase" previously aborted the
  # entire multiverse on odd-length input.
  set.seed(42)
  x <- rnorm(301)
  y <- rnorm(301)
  mv <- synchrony_multiverse(
    x, y,
    estimator = "wcc", sample_rate = 30,
    window_sec = 1, lag_sec = 0.3,
    n_surrogates = 20, surrogate_method = "phase"
  )
  expect_s3_class(mv, "bsync_multiverse")
  expect_gte(mv$robustness$n_valid, 1L)
})

test_that("generate_surrogate_circular warns and errors correctly", {
  y <- 1:15
  expect_error(
    generate_surrogate_circular(y, n_surrogates = 10, lag_max = 10),
    "too short relative to"
  )

  expect_warning(
    generate_surrogate_circular(1:30, n_surrogates = 25, lag_max = 5),
    "Limited unique shifts available"
  )
})

# =========================================================================
# --- 2. INTEGRATION TESTS (REAL PIPELINE) --------------------------------
# =========================================================================

test_that("Evaluators error on bad surrogate matrix inputs", {
  x <- 1:20
  y <- 1:20

  # Loosened the regex to ignore cli backtick formatting
  expect_error(
    wcc_surrogate(
      x,
      y,
      y_surrogates = c(1, 2, 3),
      window_size = 5,
      lag_max = 2
    ),
    "must be a matrix"
  )

  expect_error(
    wdtw_surrogate(
      x,
      y,
      y_surrogates = c(1, 2, 3),
      window_size = 5,
      lag_max = 2
    ),
    "must be a matrix"
  )

  bad_mat <- matrix(0, nrow = 10, ncol = 5)
  expect_error(
    wgranger_surrogate(x, y, y_surrogates = bad_mat, window_size = 5),
    "same number of rows"
  )
})

test_that("WCC surrogate pipeline integrates and returns valid object", {
  x <- rnorm(50)
  y <- rnorm(50)
  y_surr_mat <- generate_surrogate_circular(y, n_surrogates = 5, lag_max = 5)

  res <- wcc_surrogate(
    x,
    y,
    y_surrogates = y_surr_mat,
    window_size = 10,
    lag_max = 5
  )

  expect_s3_class(res, "wcc_surr")
  expect_type(res, "list")
  expect_equal(length(res$surrogate_z), 5)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

test_that("WDTW surrogate pipeline integrates and returns valid object", {
  x <- rnorm(50)
  y <- rnorm(50)
  y_surr_mat <- generate_surrogate_circular(y, n_surrogates = 5, lag_max = 5)

  # FIX 1: lag_max = 5 added
  res <- wdtw_surrogate(
    x,
    y,
    y_surrogates = y_surr_mat,
    window_size = 10,
    lag_max = 5
  )

  expect_s3_class(res, "wdtw_surr")
  expect_type(res, "list")
  expect_equal(length(res$surrogate_cost), 5)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

test_that("WGranger surrogate pipeline integrates and returns valid object", {
  x <- rnorm(50)
  y <- rnorm(50)
  y_surr_mat <- generate_surrogate_circular(y, n_surrogates = 5, lag_max = 5)

  # FIX 2: ar_order = 1 used instead of order = 1
  res <- wgranger_surrogate(
    x,
    y,
    y_surrogates = y_surr_mat,
    window_size = 10,
    ar_order = 1
  )

  expect_s3_class(res, "wgranger_surr")
  expect_type(res, "list")
  expect_equal(length(res$surrogate_f_xy), 5)
  expect_true(res$p_value_xy >= 0 && res$p_value_xy <= 1)
})

# =========================================================================
# --- 2b. ADDITIONAL ROBUSTNESS TESTS -------------------------------------
# =========================================================================

test_that("wcc_surrogate handles NA values in surrogate matrices without crashing", {
  set.seed(7)
  x <- rnorm(50)
  y <- rnorm(50)
  y_surr_mat <- generate_surrogate_circular(y, n_surrogates = 5, lag_max = 5)

  # Inject NAs into one surrogate column
  y_surr_mat[10:15, 2] <- NA

  # na.rm = TRUE (default): should return finite p-value
  res_true <- wcc_surrogate(
    x, y,
    y_surrogates = y_surr_mat,
    window_size = 10, lag_max = 5, na.rm = TRUE
  )
  expect_true(res_true$p_value >= 0 && res_true$p_value <= 1)

  # na.rm = FALSE: surrogate 2 will have NA-affected windows; p-value still valid
  res_false <- wcc_surrogate(
    x, y,
    y_surrogates = y_surr_mat,
    window_size = 10, lag_max = 5, na.rm = FALSE
  )
  expect_true(res_false$p_value >= 0 && res_false$p_value <= 1)
})

test_that("wcc_surrogate p-value is small for strongly coupled series", {
  set.seed(123)
  n <- 100
  x <- sin(seq(0, 4 * pi, length.out = n))
  y <- x + rnorm(n, sd = 0.05) # nearly identical to x
  y_surr_mat <- generate_surrogate_circular(y, n_surrogates = 50, lag_max = 5)

  res <- wcc_surrogate(
    x, y,
    y_surrogates = y_surr_mat,
    window_size = 20, lag_max = 5
  )

  # Observed synchrony is real; surrogate distribution should be much lower
  expect_true(res$observed_z > mean(res$surrogate_z))
  expect_true(res$p_value <= 0.10)
})

test_that("wcc_surrogate p-value is large for independent series", {
  set.seed(456)
  n <- 100
  x <- rnorm(n)
  y <- rnorm(n) # independent of x
  y_surr_mat <- generate_surrogate_circular(y, n_surrogates = 50, lag_max = 5)

  res <- wcc_surrogate(
    x, y,
    y_surrogates = y_surr_mat,
    window_size = 20, lag_max = 5
  )

  # No real synchrony; observed z should not consistently exceed surrogates
  expect_true(res$p_value > 0.05)
})

# =========================================================================
# --- 3. S3 PRINT METHOD TESTS --------------------------------------------
# =========================================================================

test_that("Print methods return silently and output text", {
  mock_wcc_obj <- list(
    observed_z = 0.8,
    surrogate_z = c(0.1, 0.2),
    p_value = 0.01,
    n_surrogates = 100,
    settings = list(statistic = "mean_abs_z")
  )
  class(mock_wcc_obj) <- c("wcc_surr", "list")

  mock_wdtw_obj <- list(
    observed_cost = 10,
    surrogate_cost = c(20, 25),
    p_value = 0.6,
    n_surrogates = 100
  )
  class(mock_wdtw_obj) <- c("wdtw_surr", "list")

  expect_invisible(print(mock_wcc_obj))

  # FIX 3: expect_message instead of capture.output
  expect_message(print(mock_wcc_obj), "significantly greater")
  expect_message(print(mock_wdtw_obj), "not significantly different")
  expect_message(print(mock_wcc_obj), "too few for stable p-values")
})

# Every "Empirical p-value: <value>" shown by a print method, read from the
# message stream (cli writes there when not interactive).
printed_p_values <- function(obj) {
  out <- paste(testthat::capture_messages(print(obj)), collapse = "")
  m <- regmatches(out, gregexpr("Empirical p-value: [^\n]+", out))[[1]]
  sub("Empirical p-value: ", "", m, fixed = TRUE)
}

mock_surr <- function(class_name, p) {
  base <- list(n_surrogates = 99999, settings = list(statistic = "mean_abs_z"))
  obj <- switch(class_name,
    wcc_surr = ,
    wphase_surr = c(base, list(
      observed_z = 0.8, surrogate_z = c(0.1, 0.2), p_value = p
    )),
    wdtw_surr = c(base, list(
      observed_cost = 10, surrogate_cost = c(20, 25), p_value = p
    )),
    wgranger_surr = c(base, list(
      observed_f_xy = 5, surrogate_f_xy = c(1, 2), p_value_xy = p[1],
      observed_f_yx = 5, surrogate_f_yx = c(1, 2), p_value_yx = p[2]
    ))
  )
  structure(obj, class = c(class_name, "list"))
}

test_that("print methods show p rounded to 4 digits, or < 0.0001", {
  # 0.123456 rounds to 0.1235. 1e-5 rounds to 0, so it shows "< 0.0001"
  # (the smallest add-one p-value with 99999 surrogates is 1 / 100000).
  # 1e-4 shows in fixed notation, not "1e-04".
  cases <- list(c(0.123456, "0.1235"), c(1e-5, "< 0.0001"), c(1e-4, "0.0001"))
  for (cls in c("wcc_surr", "wdtw_surr", "wphase_surr")) {
    for (cs in cases) {
      expect_identical(
        printed_p_values(mock_surr(cls, as.numeric(cs[1]))),
        cs[2],
        label = paste(cls, cs[1])
      )
    }
  }
  # wgranger: each direction shows its own p-value, in x -> y, y -> x order.
  expect_identical(
    printed_p_values(mock_surr("wgranger_surr", c(0.123456, 1e-5))),
    c("0.1235", "< 0.0001")
  )
  expect_identical(
    printed_p_values(mock_surr("wgranger_surr", c(1e-5, 0.123456))),
    c("< 0.0001", "0.1235")
  )
  # The display does not depend on getOption("digits").
  op <- options(digits = 3)
  on.exit(options(op), add = TRUE)
  expect_identical(printed_p_values(mock_surr("wcc_surr", 0.123456)), "0.1235")
})

test_that("a surrogate matrix with no columns aborts in every wrapper", {
  # (0 + 1) / (0 + 1) would report p = 1 as if it were a result.
  x <- rnorm(60)
  y <- rnorm(60)
  empty <- matrix(numeric(0), nrow = 60, ncol = 0)
  msg <- "must have at least one column"
  expect_error(wcc_surrogate(x, y, empty, window_size = 20, lag_max = 3), msg)
  expect_error(wdtw_surrogate(x, y, empty, window_size = 20, lag_max = 3), msg)
  expect_error(wgranger_surrogate(x, y, empty, window_size = 20), msg)
  expect_error(wphase_surrogate(x, y, empty, window_size = 20, lag_max = 3), msg)
})

# M4 acceptance-criteria tests ------------------------------------------------

test_that("M4: wcc_surrogate observed_z matches wcc() fisher_z exactly (Invariant 2)", {
  set.seed(1)
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B
  y_surr <- generate_surrogate_circular(y, n_surrogates = 5, lag_max = 10)

  # mean_abs_z path
  res_wcc_maz <- wcc(x, y, window_size = 96, lag_max = 10, statistic = "mean_abs_z")
  res_surr_maz <- wcc_surrogate(
    x, y,
    y_surrogates = y_surr,
    window_size = 96, lag_max = 10,
    statistic = "mean_abs_z"
  )
  expect_equal(res_surr_maz$observed_z, res_wcc_maz$aggregate[[1]])

  # peak path
  res_wcc_peak <- wcc(x, y, window_size = 96, lag_max = 10, statistic = "peak")
  res_surr_peak <- wcc_surrogate(
    x, y,
    y_surrogates = y_surr,
    window_size = 96, lag_max = 10,
    statistic = "peak"
  )
  expect_equal(res_surr_peak$observed_z, res_wcc_peak$aggregate[[1]])
})

test_that("M4: surrogate loop aggregate equals observed path (cross-path Invariant 2)", {
  # Pass y itself as the sole surrogate; surrogate_z[1] must equal observed_z.
  # This exercises the surrogate loop's wcc_aggregate(grid_df$row) call site
  # independently from the observed wcc_aggregate(results_df$i) call site —
  # the two code paths use different group labels but the same partition.
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B
  y_self <- matrix(y, ncol = 1)

  for (stat in c("mean_abs_z", "peak")) {
    res_obs <- wcc(x, y, window_size = 96, lag_max = 10, statistic = stat)
    res_surr <- wcc_surrogate(
      x, y,
      y_surrogates = y_self,
      window_size = 96, lag_max = 10,
      statistic = stat
    )
    expect_equal(
      res_surr$surrogate_z[[1]], res_obs$aggregate[[1]],
      label = paste0("surrogate loop == observed for statistic='", stat, "'")
    )
  }
})

test_that("M4: wcc_surrogate statistic arg is validated and recorded", {
  set.seed(2)
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B
  y_surr <- generate_surrogate_circular(y, n_surrogates = 3, lag_max = 10)

  res_peak <- wcc_surrogate(
    x, y,
    y_surrogates = y_surr,
    window_size = 96, lag_max = 10,
    statistic = "peak"
  )
  expect_equal(res_peak$settings$statistic, "peak")

  expect_error(
    wcc_surrogate(
      x, y,
      y_surrogates = y_surr,
      window_size = 96, lag_max = 10,
      statistic = "bad"
    ),
    "should be one of"
  )
})

test_that("M4: peak p-value is small for coupled series, large for independent", {
  set.seed(123)
  n <- 100
  x_coupled <- sin(seq(0, 4 * pi, length.out = n))
  y_coupled <- x_coupled + rnorm(n, sd = 0.05)
  y_surr_coupled <- generate_surrogate_circular(y_coupled,
    n_surrogates = 50,
    lag_max = 5
  )

  res_coupled <- wcc_surrogate(
    x_coupled, y_coupled,
    y_surrogates = y_surr_coupled,
    window_size = 20, lag_max = 5,
    statistic = "peak"
  )
  expect_true(res_coupled$observed_z > mean(res_coupled$surrogate_z))
  expect_true(res_coupled$p_value <= 0.10)

  set.seed(456)
  x_indep <- rnorm(n)
  y_indep <- rnorm(n)
  y_surr_indep <- generate_surrogate_circular(y_indep,
    n_surrogates = 50,
    lag_max = 5
  )

  res_indep <- wcc_surrogate(
    x_indep, y_indep,
    y_surrogates = y_surr_indep,
    window_size = 20, lag_max = 5,
    statistic = "peak"
  )
  expect_true(res_indep$p_value > 0.05)
})


# M5 AC4: surrogate engine tests ----------------------------------------------

test_that("AC4: run_surrogate_engine returns correct shape", {
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B
  y_surr <- generate_surrogate_circular(y, n_surrogates = 3, lag_max = 10)

  grid <- build_surface_grid(
    n_x = length(x), window_size = 96, window_increment = 1,
    lag_max = 10, lag_increment = 1, lagged = TRUE
  )
  x_cpp <- as.double(x)

  # Scalar aggregate (WCC-like)
  compute_scalar <- function(xv, y_col, g) {
    mean(abs(r_to_z(calc_wcc_cpp(
      x = xv, y = as.double(y_col),
      i_vals = g$i_vals, tau_vals = g$tau_vals,
      w_max = g$w_max, na_rm = TRUE
    ))), na.rm = TRUE)
  }

  result <- run_surrogate_engine(x_cpp, y_surr, grid, compute_scalar, numeric(1))
  expect_true(is.numeric(result))
  expect_equal(length(result), 3L)
  expect_false(is.data.frame(result)) # aggregate-only: no results_df
})

test_that("AC4: surrogate result objects carry no results_df (Invariant 7)", {
  # The aggregate-only path must never materialize a per-cell results_df, even
  # in the returned object. Assert all three surrogate objects omit it.
  # A 600-sample slice exercises the same code path as the full series — this
  # is a structural assertion, not a numeric one — at a fraction of the WDTW
  # cost (DTW is O(window^2) per cell).
  x <- sim_dyad$x_A[1:600]
  y <- sim_dyad$x_B[1:600]
  y_surr <- generate_surrogate_circular(y, n_surrogates = 3, lag_max = 10)

  res_wcc <- wcc_surrogate(x, y, y_surrogates = y_surr, window_size = 96, lag_max = 10)
  res_wdtw <- wdtw_surrogate(x, y, y_surrogates = y_surr, window_size = 96, lag_max = 10)
  res_wg <- wgranger_surrogate(x, y, y_surrogates = y_surr, window_size = 96)

  expect_null(res_wcc$results_df)
  expect_null(res_wdtw$results_df)
  expect_null(res_wg$results_df)
})

test_that("AC4: seeded surrogate p-values are reproducible (regression guard)", {
  # Frozen against the current implementation on a fixed seed; guards against a
  # tail-direction flip or off-by-one in the empirical p-value count. AR(1)
  # series give a non-boundary p-value that actually exercises tail counting.
  # Add-one form (b + 1) / (n + 1), Phipson & Smyth (2010, p. 6); the b counts
  # (91 and 83 of 99) are unchanged from the earlier b / n pin.
  set.seed(20260629)
  n <- 120
  x <- as.numeric(stats::arima.sim(list(ar = 0.5), n))
  y <- as.numeric(stats::arima.sim(list(ar = 0.5), n))

  y_surr <- generate_surrogate_circular(y, n_surrogates = 99, lag_max = 5)
  res_wcc <- wcc_surrogate(x, y, y_surrogates = y_surr, window_size = 30, lag_max = 5)
  expect_equal(res_wcc$p_value, (91 + 1) / (99 + 1), tolerance = 1e-12)

  res_wdtw <- wdtw_surrogate(x, y, y_surrogates = y_surr, window_size = 30, lag_max = 5)
  expect_equal(res_wdtw$p_value, (83 + 1) / (99 + 1), tolerance = 1e-12)
})

# Add-one p-value, Phipson & Smyth (2010, p. 6): p = (b + 1) / (n + 1), where
# n is the number of surrogate columns and b counts surrogate statistics at
# least as extreme as the observed one. b is counted here with a plain loop,
# independent of the wrappers' vectorized sum.
count_at_least_as_extreme <- function(observed, surrogates, tail) {
  b <- 0
  for (s in surrogates) {
    if (tail == "upper" && s >= observed) b <- b + 1
    if (tail == "lower" && s <= observed) b <- b + 1
  }
  b
}

# Run all four wrappers on one dyad and 19 circular-shift surrogates; return
# each statistic's observed value, surrogate draws, p-value, and tail.
add_one_probe <- function(x, y, seed) {
  set.seed(seed)
  ys <- generate_surrogate_circular(y, n_surrogates = 19, lag_max = 4)
  args <- list(window_size = 40, lag_max = 4, window_increment = 10)
  wc <- do.call(wcc_surrogate, c(list(x, y, ys), args))
  wd <- do.call(wdtw_surrogate, c(list(x, y, ys), args))
  wg <- wgranger_surrogate(x, y, ys, window_size = 40, window_increment = 10)
  wp <- do.call(wphase_surrogate, c(list(x, y, ys), args))
  list(
    n = ncol(ys),
    stats = list(
      wcc = list(wc$observed_z, wc$surrogate_z, wc$p_value, "upper"),
      wdtw = list(wd$observed_cost, wd$surrogate_cost, wd$p_value, "lower"),
      wgranger_xy = list(wg$observed_f_xy, wg$surrogate_f_xy, wg$p_value_xy, "upper"),
      wgranger_yx = list(wg$observed_f_yx, wg$surrogate_f_yx, wg$p_value_yx, "upper"),
      wphase = list(wp$observed_z, wp$surrogate_z, wp$p_value, "upper")
    )
  )
}

check_add_one <- function(probe, b_ok) {
  for (nm in names(probe$stats)) {
    s <- probe$stats[[nm]]
    b <- count_at_least_as_extreme(s[[1]], s[[2]], s[[4]])
    expect_true(b_ok(b), label = paste(nm, "b =", b))
    expect_equal(s[[3]], (b + 1) / (probe$n + 1), tolerance = 1e-12, label = nm)
  }
}

test_that("all four wrappers return (b + 1) / (n + 1) when b = 0", {
  # Bidirectional VAR(1) coupling (x <-> y, 0.6 each way): every statistic,
  # both Granger directions included, beats all 19 shifted surrogates.
  set.seed(1)
  n <- 200
  e1 <- rnorm(n)
  e2 <- rnorm(n)
  x <- numeric(n)
  y <- numeric(n)
  for (t in 2:n) {
    x[t] <- 0.6 * y[t - 1] + e1[t]
    y[t] <- 0.6 * x[t - 1] + e2[t]
  }
  probe <- add_one_probe(x, y, seed = 2)
  check_add_one(probe, function(b) b == 0)
  # b = 0 gives the smallest possible p-value, 1 / 20, never 0.
  expect_equal(probe$stats$wcc[[3]], 1 / 20, tolerance = 1e-12)
})

test_that("all four wrappers return (b + 1) / (n + 1) when 0 < b < n", {
  set.seed(2)
  x <- rnorm(200)
  y <- rnorm(200)
  probe <- add_one_probe(x, y, seed = 102)
  check_add_one(probe, function(b) b > 0 && b < 19)
})

test_that("AC4: WDTW fast_method evaluates observed windows at lag 0 (no over-count)", {
  # Regression guard for the fast-path grid: surrogates must cover exactly the
  # observed lagged surface's window positions, evaluated at tau = 0 — not a
  # lag-free grid that shifts windows past the series end and over-counts.
  # The window-alignment identity is length-independent; a 600-sample slice
  # keeps the guard while cutting the WDTW cost.
  x <- sim_dyad$x_A[1:600]
  y <- sim_dyad$x_B[1:600]
  y_self <- matrix(y, ncol = 1)

  obs <- wdtw(x, y, window_size = 96, lag_max = 10, scale_method = "none")
  obs_lag0_mean <- mean(
    obs$results_df$dtw_dist[obs$results_df$tau == 0],
    na.rm = TRUE
  )

  res_fast <- wdtw_surrogate(
    x, y,
    y_surrogates = y_self,
    window_size = 96, lag_max = 10, scale_method = "none",
    fast_method = TRUE
  )

  # y is its own surrogate ⇒ fast cost == observed mean over the same windows
  # at lag 0. If the fast grid over-counted (extra out-of-range windows), this
  # would diverge.
  expect_equal(res_fast$surrogate_cost[[1]], obs_lag0_mean, tolerance = 1e-12)
})

test_that("AC4: run_surrogate_engine supports named-numeric return (Granger-like)", {
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B
  y_surr <- generate_surrogate_circular(y, n_surrogates = 3, lag_max = 10)

  grid <- build_surface_grid(
    n_x = length(x), window_size = 96, window_increment = 1,
    lagged = FALSE
  )
  x_cpp <- as.double(x)

  compute_pair <- function(xv, y_col, g) {
    s <- calc_wgranger_cpp(xv, as.double(y_col), g$i_vals, g$w_max, p = 1L)
    c(f_xy = mean(s$f_xy, na.rm = TRUE), f_yx = mean(s$f_yx, na.rm = TRUE))
  }

  result <- run_surrogate_engine(x_cpp, y_surr, grid, compute_pair, numeric(2))
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2L, 3L))
  expect_equal(rownames(result), c("f_xy", "f_yx"))
})

test_that("AC4: WDTW cross-path Invariant-2 (y as its sole surrogate)", {
  # Pass y itself as the sole surrogate at lag=0 path (scale_method='none').
  # The surrogate's mean DTW cost must equal wdtw()$aggregate[['mean_distance']].
  # This is an exact algebraic identity, independent of series length, so a
  # 600-sample slice preserves it at a fraction of the WDTW cost.
  x <- sim_dyad$x_A[1:600]
  y <- sim_dyad$x_B[1:600]
  y_self <- matrix(y, ncol = 1)

  res_obs <- wdtw(x, y, window_size = 96, lag_max = 10, scale_method = "none")
  res_surr <- wdtw_surrogate(
    x, y,
    y_surrogates = y_self,
    window_size = 96, lag_max = 10, scale_method = "none"
  )

  expect_equal(
    res_surr$surrogate_cost[[1]],
    res_obs$aggregate[["mean_distance"]],
    tolerance = 1e-12
  )
})

test_that("AC4: Granger cross-path Invariant-2 (y as its sole surrogate)", {
  # Pass y itself as the sole surrogate.
  # surrogate_f_xy[1] must equal observed_f_xy, and same for f_yx.
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B
  y_self <- matrix(y, ncol = 1)

  res_obs <- wgranger(x, y, window_size = 96)
  res_surr <- wgranger_surrogate(
    x, y,
    y_surrogates = y_self,
    window_size = 96
  )

  expect_equal(res_surr$surrogate_f_xy[[1]], res_obs$aggregate[["f_xy"]], tolerance = 1e-12)
  expect_equal(res_surr$surrogate_f_yx[[1]], res_obs$aggregate[["f_yx"]], tolerance = 1e-12)
})

# M7 Phase B regression tests -------------------------------------------------

test_that("M7: print.wcc_surr labels have no double-colon (cli regression)", {
  set.seed(42)
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B
  y_surr <- generate_surrogate_circular(y, n_surrogates = 20, lag_max = 10)

  res_maz <- wcc_surrogate(x, y,
    y_surrogates = y_surr, window_size = 96, lag_max = 10,
    statistic = "mean_abs_z"
  )
  res_peak <- wcc_surrogate(x, y,
    y_surrogates = y_surr, window_size = 96, lag_max = 10,
    statistic = "peak"
  )

  # Labels must appear exactly once, with a single colon (no ": :" double colon)
  expect_message(print(res_maz), "Observed Mean Abs\\. Fisher")
  expect_message(print(res_maz), "Average Null Mean Abs\\. Fisher")
  expect_message(print(res_peak), "Observed Mean Peak Abs\\. Fisher")
  expect_message(print(res_peak), "Average Null Mean Peak Abs\\. Fisher")

  # Capture full output and assert absence of double colon
  out_maz <- capture.output(suppressMessages(
    withCallingHandlers(print(res_maz), message = function(m) {
      cat(conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  ))
  expect_false(any(grepl(": :", out_maz)),
    label = "print.wcc_surr (mean_abs_z) must not contain double colon ': :'"
  )

  out_peak <- capture.output(suppressMessages(
    withCallingHandlers(print(res_peak), message = function(m) {
      cat(conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  ))
  expect_false(any(grepl(": :", out_peak)),
    label = "print.wcc_surr (peak) must not contain double colon ': :'"
  )
})


# wphase surrogate: Invariant 2, Invariant 6, and calibration -----------------

test_that("Invariant 2: wphase_surrogate observed_z equals wphase aggregate exactly", {
  set.seed(11)
  x <- sim_dyad$z_A[1:300]
  y <- sim_dyad$z_B[1:300]
  ys <- generate_surrogate_circular(y, n_surrogates = 5)

  obs <- wphase(x, y, window_size = 64, lag_max = 5, window_increment = 8)
  surr <- wphase_surrogate(x, y,
    y_surrogates = ys,
    window_size = 64, lag_max = 5, window_increment = 8
  )

  expect_identical(surr$observed_z, obs$aggregate[[1]])

  # Cross-path check: y itself as the sole "surrogate" must reproduce the
  # observed aggregate through the aggregate-only surrogate path.
  self <- wphase_surrogate(x, y,
    y_surrogates = matrix(y, ncol = 1),
    window_size = 64, lag_max = 5, window_increment = 8
  )
  expect_equal(self$surrogate_z[1], obs$aggregate[[1]], tolerance = 1e-12)
})

test_that("Invariant 6: wphase_surrogate is reproducible under set.seed", {
  x <- sim_dyad$z_A[1:300]
  y <- sim_dyad$z_B[1:300]

  run_once <- function() {
    set.seed(99)
    ys <- generate_surrogate_circular(y, n_surrogates = 19)
    wphase_surrogate(x, y,
      y_surrogates = ys,
      window_size = 64, lag_max = 5, window_increment = 8
    )
  }
  a <- run_once()
  b <- run_once()
  expect_identical(a$p_value, b$p_value)
  expect_identical(a$surrogate_z, b$surrogate_z)
})

test_that("wphase surrogate p-values are calibrated under the null (Type I)", {
  skip_on_cran()
  # 200 seeded replicates of independent white-noise pairs, n_surrogates = 99.
  # With B = 99 surrogates, p = (b + 1) / 100 <= .05 iff at most 4 surrogates
  # >= observed, which has probability exactly 5/100 under exchangeability. Over R = 200
  # replicates the rejection count is Binomial(200, .05): mean 10,
  # SE = sqrt(200 * .05 * .95) = 3.08; a +/- 3 SE band is 10 +/- 9.25,
  # so integer counts in [1, 19].
  set.seed(20260829)
  rep_seeds <- sample.int(.Machine$integer.max, 200)
  pvals <- vapply(rep_seeds, function(s) {
    set.seed(s)
    x <- rnorm(300)
    y <- rnorm(300)
    ys <- generate_surrogate_circular(y, n_surrogates = 99)
    wphase_surrogate(x, y,
      y_surrogates = ys,
      window_size = 64, lag_max = 5, window_increment = 8
    )$p_value
  }, numeric(1))
  rejections <- sum(pvals <= 0.05)
  expect_gte(rejections, 1)
  expect_lte(rejections, 19)
})

test_that("wphase surrogate test has power against phase-locked signals", {
  skip_on_cran()
  # Coupling design: a shared frequency-wandering sinusoid (theta = 8 Hz
  # carrier at fs = 128 plus a cumulative-normal phase wander, sd = 0.15/step)
  # with independent additive noise, sd = 0.5. PLV is invariant to constant
  # phase offsets, so a strictly periodic common component would give the
  # circular-shift null NO power (the shift only offsets the phase); the
  # shared wander is what a shift misaligns. At this design a 60-replicate
  # pilot rejected at rate .983 (~.95+ expected power); the assertion floor
  # is .80, leaving ~6 SE of Monte-Carlo headroom at power .95
  # (SE = sqrt(.95*.05/200) = .0154).
  set.seed(20260830)
  rep_seeds <- sample.int(.Machine$integer.max, 200)
  pvals <- vapply(rep_seeds, function(s) {
    set.seed(s)
    n <- 300
    t <- 0:(n - 1)
    theta <- 2 * pi * 8 * t / 128 + cumsum(rnorm(n, sd = 0.15))
    x <- sin(theta) + rnorm(n, sd = 0.5)
    y <- sin(theta + 0.8) + rnorm(n, sd = 0.5)
    ys <- generate_surrogate_circular(y, n_surrogates = 99)
    wphase_surrogate(x, y,
      y_surrogates = ys,
      window_size = 64, lag_max = 5, window_increment = 8
    )$p_value
  }, numeric(1))
  expect_gte(mean(pvals <= 0.05), 0.80)
})

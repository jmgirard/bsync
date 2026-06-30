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
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B
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
  set.seed(20260629)
  n <- 120
  x <- as.numeric(stats::arima.sim(list(ar = 0.5), n))
  y <- as.numeric(stats::arima.sim(list(ar = 0.5), n))

  y_surr <- generate_surrogate_circular(y, n_surrogates = 99, lag_max = 5)
  res_wcc <- wcc_surrogate(x, y, y_surrogates = y_surr, window_size = 30, lag_max = 5)
  expect_equal(res_wcc$p_value, 91 / 99, tolerance = 1e-12)

  res_wdtw <- wdtw_surrogate(x, y, y_surrogates = y_surr, window_size = 30, lag_max = 5)
  expect_equal(res_wdtw$p_value, 83 / 99, tolerance = 1e-12)
})

test_that("AC4: WDTW fast_method evaluates observed windows at lag 0 (no over-count)", {
  # Regression guard for the fast-path grid: surrogates must cover exactly the
  # observed lagged surface's window positions, evaluated at tau = 0 — not a
  # lag-free grid that shifts windows past the series end and over-counts.
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B
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
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B
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

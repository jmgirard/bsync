# AC0 preflight: external-oracle validation (DESIGN.md §13 Layer 2)
#
# These tests freeze golden values produced by external reference packages
# (dtw v1.23-3, lmtest v0.9-40) against tiny fixed inputs. The assertions
# run against the *committed frozen numbers*, not against the live packages,
# so there is no CI dependency on dtw/lmtest.
#
# Provenance:
#   dtw golden (WDTW L1): computed 2026-06-29 with dtw::dtw(x[1:10], y[1:10],
#     step.pattern = dtw::symmetric1, distance.only = TRUE)$distance
#     where x = sin(seq(0,2*pi,20)), y = cos(seq(0,2*pi,20)), set.seed(0)
#   Pure-R oracle (WDTW L2): dtw_oracle_l2() below — squared-difference DP,
#     confirmed bsync matches to <1e-12.
#   lmtest golden (Granger): computed 2026-06-29 with
#     lmtest::grangertest(y ~ x, order=1) and lmtest::grangertest(x ~ y, order=1)
#     on 30-point series (set.seed(42)), treating window_size = length(x) = 30
#     so bsync produces a single window matching lmtest's full-series test.

# Frozen fixture inputs -------------------------------------------------------

oracle_x10 <- sin(seq(0, 2 * pi, length.out = 20))[1:10]
oracle_y10 <- cos(seq(0, 2 * pi, length.out = 20))[1:10]

granger_series <- local({
  set.seed(42)
  n <- 30
  x <- rnorm(n)
  list(
    x = x,
    y = 0.8 * dplyr::lag(x, 1, default = 0) + rnorm(n, sd = 0.3)
  )
})


# WDTW Layer-2 oracle tests ---------------------------------------------------

test_that("AC0: WDTW L1 matches dtw::dtw (symmetric1) golden value", {
  # dtw golden: 6.221439863967936  (dtw v1.23-3, symmetric1 step pattern)
  dtw_golden_l1 <- 6.221439863967936

  result <- calc_wdtw_cpp(
    x = as.double(oracle_x10),
    y = as.double(oracle_y10),
    i_vals = 1L,
    tau_vals = 0L,
    w_max = 9L,
    use_l2 = FALSE,
    local_scale = FALSE
  )

  expect_equal(result[[1]], dtw_golden_l1, tolerance = 1e-10)
})

test_that("AC0: WDTW L2 matches pure-R squared-difference DTW oracle", {
  # Pure-R oracle: symmetric1 DP with local cost (x_i - y_j)^2
  dtw_oracle_l2 <- function(x, y) {
    n <- length(x)
    D <- matrix(Inf, nrow = n + 1, ncol = n + 1)
    D[1, 1] <- 0
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        D[i + 1, j + 1] <- (x[i] - y[j])^2 +
          min(D[i, j], D[i + 1, j], D[i, j + 1])
      }
    }
    D[n + 1, n + 1]
  }

  oracle_val <- dtw_oracle_l2(oracle_x10, oracle_y10)

  result <- calc_wdtw_cpp(
    x = as.double(oracle_x10),
    y = as.double(oracle_y10),
    i_vals = 1L,
    tau_vals = 0L,
    w_max = 9L,
    use_l2 = TRUE,
    local_scale = FALSE
  )

  expect_equal(result[[1]], oracle_val, tolerance = 1e-12)
})


# Granger Layer-2 oracle tests ------------------------------------------------

test_that("AC0: Granger F-stats match lmtest::grangertest golden values", {
  # lmtest golden (lmtest v0.9-40, grangertest(y ~ x, order=1)):
  #   f_xy = 287.6196011323,  p_xy = 1.409584e-15
  #   f_yx =   0.2771169716,  p_yx = 6.030594e-01
  golden_f_xy <- 287.6196011323
  golden_p_xy <- 1.409584053095196e-15
  golden_f_yx <- 0.2771169716
  golden_p_yx <- 6.030594127984489e-01

  res <- wgranger(
    x = granger_series$x,
    y = granger_series$y,
    window_size = 30,
    ar_order = 1
  )
  df <- res$results_df

  expect_equal(df$f_xy[1], golden_f_xy, tolerance = 1e-6)
  expect_equal(df$p_xy[1], golden_p_xy, tolerance = 1e-10)
  expect_equal(df$f_yx[1], golden_f_yx, tolerance = 1e-6)
  expect_equal(df$p_yx[1], golden_p_yx, tolerance = 1e-6)
})

test_that("AC0: Granger R-level oracle (lm-based) matches C++ core", {
  # Independent pure-R implementation using lm() — same AR structure as the C++ core.
  # This is the pure-R oracle analogous to wcc_oracle(), decoupled from lmtest.
  granger_oracle <- function(x, y, window_size, ar_order = 1) {
    p <- ar_order
    n_eff <- window_size - p

    y_tgt <- y[(p + 1):window_size]
    x_tgt <- x[(p + 1):window_size]

    # Build lag matrices
    Y_lags <- sapply(seq_len(p), function(k) y[(p + 1 - k):(window_size - k)])
    X_lags <- sapply(seq_len(p), function(k) x[(p + 1 - k):(window_size - k)])

    df_u <- as.data.frame(cbind(y_tgt = y_tgt, Y_lags, X_lags))
    df_r_y <- as.data.frame(cbind(y_tgt = y_tgt, Y_lags))
    df_r_x <- as.data.frame(cbind(x_tgt = x_tgt, X_lags))
    df_u_x <- as.data.frame(cbind(x_tgt = x_tgt, Y_lags, X_lags))

    f_xy <- {
      rss_u <- sum(resid(lm(y_tgt ~ ., data = df_u))^2)
      rss_r <- sum(resid(lm(y_tgt ~ ., data = df_r_y))^2)
      df2 <- n_eff - 2 * p - 1
      ((rss_r - rss_u) / p) / (rss_u / df2)
    }
    f_yx <- {
      rss_u <- sum(resid(lm(x_tgt ~ ., data = df_u_x))^2)
      rss_r <- sum(resid(lm(x_tgt ~ ., data = df_r_x))^2)
      df2 <- n_eff - 2 * p - 1
      ((rss_r - rss_u) / p) / (rss_u / df2)
    }

    list(f_xy = f_xy, f_yx = f_yx)
  }

  oracle <- granger_oracle(
    x = granger_series$x,
    y = granger_series$y,
    window_size = 30,
    ar_order = 1
  )

  res <- wgranger(
    x = granger_series$x,
    y = granger_series$y,
    window_size = 30,
    ar_order = 1
  )

  expect_equal(res$results_df$f_xy[1], oracle$f_xy, tolerance = 1e-9)
  expect_equal(res$results_df$f_yx[1], oracle$f_yx, tolerance = 1e-9)
})

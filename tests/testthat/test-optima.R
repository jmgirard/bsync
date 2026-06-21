# Helper function to quickly mock up a wcc_res object for testing
create_mock_wcc <- function(wcc_vals, tau_max, i_val = 100) {
  lags <- seq(-tau_max, tau_max)
  df <- data.frame(
    i = i_val,
    tau = lags,
    wcc = wcc_vals
  )

  obj <- list(
    results_df = df,
    settings = list(lag_max = tau_max)
  )
  class(obj) <- c("wcc_res", "list")
  obj
}

test_that("pick_optima validates inputs and handles errors properly", {
  mock_wcc <- create_mock_wcc(c(0.1, 0.5, 0.9, 0.5, 0.1), tau_max = 2)

  # Test class check
  expect_error(
    pick_optima(list(), search_method = "global"),
    "<wcc_res> or <wdtw_res>"
  )

  # Test even L_size C++ exception (requires local search)
  expect_error(
    pick_optima(mock_wcc, L_size = 4, search_method = "local"),
    "L_size must be an odd integer"
  )

  # Test missing L_size for local search
  expect_error(
    pick_optima(mock_wcc, search_method = "local"),
    "must be provided"
  )
})

test_that("pick_optima finds a perfect central peak (local search)", {
  # Peak exactly at lag 0
  mock_wcc <- create_mock_wcc(c(0.1, 0.5, 0.9, 0.5, 0.1), tau_max = 2)
  res <- pick_optima(mock_wcc, L_size = 3, search_method = "local")

  expect_equal(res$optimum_lag, 0)
  expect_equal(res$optimum_value, 0.9)
  expect_equal(res$i, 100)
})

test_that("pick_optima finds shifted off-center peaks (local search)", {
  # Peak at lag +1
  mock_pos <- create_mock_wcc(c(0.1, 0.2, 0.5, 0.8, 0.3), tau_max = 2)
  res_pos <- pick_optima(mock_pos, L_size = 3, search_method = "local")

  expect_equal(res_pos$optimum_lag, 1)
  expect_equal(res_pos$optimum_value, 0.8)

  # Peak at lag -2 (requires tau_max >= 3 so lag -3 exists for the flank)
  mock_neg <- create_mock_wcc(c(0.4, 0.9, 0.5, 0.2, 0.1, 0.1, 0.1), tau_max = 3)
  res_neg <- pick_optima(mock_neg, L_size = 3, search_method = "local")

  expect_equal(res_neg$optimum_lag, -2)
  expect_equal(res_neg$optimum_value, 0.9)
})

test_that("strict_monotonic toggle enforces flank constraints (local search)", {
  # Peak is at lag 0 (0.9).
  # Moving left from the peak, the values go from 0.5 to 0.7.
  # This violates the rule that values must strictly decrease away from the peak.
  mock_wcc <- create_mock_wcc(c(0.7, 0.5, 0.9, 0.4, 0.1), tau_max = 2)

  # With strict_monotonic = FALSE (default), it should find the local max
  res_loose <- pick_optima(
    mock_wcc,
    L_size = 5,
    strict_monotonic = FALSE,
    search_method = "local"
  )
  expect_equal(res_loose$optimum_lag, 0)

  # With strict_monotonic = TRUE, it should reject this peak and return NAs
  res_strict <- pick_optima(
    mock_wcc,
    L_size = 5,
    strict_monotonic = TRUE,
    search_method = "local"
  )
  expect_true(is.na(res_strict$optimum_lag))
  expect_true(is.na(res_strict$optimum_value))
})

test_that("pick_optima finds the correct absolute extremum (global search)", {
  # Create a WCC surface with two peaks, one higher than the other
  mock_wcc <- create_mock_wcc(c(0.2, 0.9, 0.4, 0.7, 0.1), tau_max = 2)

  # Global search should bypass L_size/flanks and just grab the absolute max
  res_global <- pick_optima(mock_wcc, search_method = "global")

  expect_equal(res_global$optimum_lag, -1)
  expect_equal(res_global$optimum_value, 0.9)
})

# Helper function to quickly mock up a wdtw_res object for testing
create_mock_wdtw <- function(dist_vals, tau_max, i_val = 100) {
  lags <- seq(-tau_max, tau_max)
  df <- data.frame(
    i = i_val,
    tau = lags,
    dtw_dist = dist_vals
  )

  obj <- list(
    results_df = df,
    settings = list(lag_max = tau_max)
  )
  class(obj) <- c("wdtw_res", "list")
  obj
}

test_that("print.wcc_optima produces expected console output", {
  mock_wcc <- create_mock_wcc(c(0.1, 0.5, 0.9, 0.5, 0.1), tau_max = 2)
  res <- pick_optima(mock_wcc, L_size = 3, search_method = "local")

  # This will capture the full cli header and the dataframe
  expect_snapshot(print(res))
})

test_that("summary.wcc_optima calculates and prints correctly", {
  # Create a scenario with one valid peak and one invalid peak
  df_valid <- data.frame(i = 100, tau = -2:2, wcc = c(0.1, 0.5, 0.9, 0.5, 0.1))
  df_invalid <- data.frame(
    i = 101,
    tau = -2:2,
    wcc = c(0.9, 0.7, 0.5, 0.4, 0.1)
  )

  mock_wcc <- list(
    results_df = rbind(df_valid, df_invalid),
    settings = list(lag_max = 2)
  )
  class(mock_wcc) <- c("wcc_res", "list")

  res <- pick_optima(
    mock_wcc,
    L_size = 5,
    strict_monotonic = TRUE,
    search_method = "local"
  )

  expect_snapshot(summary(res))
})

test_that("print.wdtw_optima produces expected console output", {
  mock_wdtw <- create_mock_wdtw(c(5.0, 3.2, 1.1, 3.2, 5.0), tau_max = 2)
  res <- pick_optima(mock_wdtw)

  expect_snapshot(print(res))
})

test_that("summary.wdtw_optima calculates and prints correctly", {
  # Mock a negative lag (y leads x) scenario
  mock_wdtw <- create_mock_wdtw(c(1.1, 3.2, 5.0, 6.1, 7.0), tau_max = 2)
  res <- pick_optima(mock_wdtw, search_method = "global")

  expect_snapshot(summary(res))
})

test_that("wcc accurately calculates correlations for identical and shifted series", {
  # Use a continuous sine wave to prevent zero-variance edge cases and boundary limits
  master <- sin(seq(0, 10, length.out = 50))

  # 1. Identical series test
  x_ident <- master[1:30]
  y_ident <- master[1:30]

  res_ident <- wcc(x_ident, y_ident, window_size = 10, lag_max = 3)
  df_ident <- res_ident$results_df

  # The correlation at lag 0 should be precisely 1 for all time windows
  lag_zero_wcc <- df_ident$wcc[df_ident$tau == 0]
  expect_equal(lag_zero_wcc, rep(1, length(lag_zero_wcc)))

  # 2. Shifted series test (y is delayed by 1 step relative to x)
  x_shift <- master[2:31]
  y_shift <- master[1:30]

  res_shift <- wcc(x_shift, y_shift, window_size = 10, lag_max = 3)
  opt_shift <- pick_optima(res_shift, search_method = "global")

  # The global peak for every window should be exactly 1 at a lag of magnitude 1.
  expect_equal(abs(opt_shift$optimum_lag), rep(1, nrow(opt_shift)))
  expect_equal(opt_shift$optimum_value, rep(1, nrow(opt_shift)))
})

test_that("wdtw accurately calculates distances for identical and shifted series", {
  master <- sin(seq(0, 10, length.out = 50))

  # 1. Identical series test
  x_ident <- master[1:30]
  y_ident <- master[1:30]

  # Set scale_method = "none" to verify the absolute raw distance logic
  res_ident <- wdtw(
    x_ident,
    y_ident,
    window_size = 10,
    lag_max = 3,
    scale_method = "none"
  )
  df_ident <- res_ident$results_df

  # The distance at lag 0 should be exactly 0 for all time windows
  lag_zero_dtw <- df_ident$dtw_dist[df_ident$tau == 0]
  expect_equal(lag_zero_dtw, rep(0, length(lag_zero_dtw)))

  # 2. Shifted series test (y is delayed by 1 step relative to x)
  x_shift <- master[2:31]
  y_shift <- master[1:30]

  res_shift <- wdtw(
    x_shift,
    y_shift,
    window_size = 10,
    lag_max = 3,
    scale_method = "none"
  )
  opt_shift <- pick_optima(res_shift, search_method = "global")

  # The global minimum distance should be exactly 0 at a lag of magnitude 1
  expect_equal(abs(opt_shift$optimum_lag), rep(1, nrow(opt_shift)))
  expect_equal(opt_shift$optimum_value, rep(0, nrow(opt_shift)))
})

test_that("pick_optima rejects invalid search_method arguments", {
  mock_wcc <- create_mock_wcc(c(0.1, 0.5, 0.9, 0.5, 0.1), tau_max = 2)

  expect_error(
    pick_optima(mock_wcc, search_method = "magic"),
    "must be either"
  )
})

test_that("pick_optima correctly applies thresholds to filter weak optima", {
  # 1. WCC Threshold (Optima below threshold become NA)
  # The peak is 0.3. A threshold of 0.5 should wipe it out.
  mock_wcc <- create_mock_wcc(c(0.1, 0.2, 0.3, 0.2, 0.1), tau_max = 2)
  res_wcc <- pick_optima(mock_wcc, search_method = "global", threshold = 0.5)

  expect_true(is.na(res_wcc$optimum_lag))
  expect_true(is.na(res_wcc$optimum_value))

  # 2. WDTW Threshold (Optima distance above threshold become NA)
  # The minimum distance is 5.0. A threshold of 3.0 should wipe it out.
  mock_wdtw <- create_mock_wdtw(c(10.0, 8.0, 5.0, 8.0, 10.0), tau_max = 2)
  res_wdtw <- pick_optima(mock_wdtw, search_method = "global", threshold = 3.0)

  expect_true(is.na(res_wdtw$optimum_lag))
  expect_true(is.na(res_wdtw$optimum_value))
})

test_that("print.wcc_optima handles empty objects and row limits correctly", {
  # Test with > 5 rows to trigger the "remaining rows" text block
  df_large <- data.frame(i = 1:10, tau = rep(0, 10), wcc = rep(0.9, 10))
  mock_wcc_large <- list(results_df = df_large, settings = list(lag_max = 2))
  class(mock_wcc_large) <- c("wcc_res", "list")

  res_large <- pick_optima(mock_wcc_large, search_method = "global")
  expect_snapshot(print(res_large))

  # Test with an empty dataframe (0 rows)
  df_empty <- data.frame(i = integer(), tau = integer(), wcc = numeric())
  mock_wcc_empty <- list(results_df = df_empty, settings = list(lag_max = 2))
  class(mock_wcc_empty) <- c("wcc_res", "list")

  res_empty <- pick_optima(mock_wcc_empty, search_method = "global")
  expect_snapshot(print(res_empty))
})

test_that("summary methods gracefully handle objects with zero valid optima", {
  mock_wcc <- create_mock_wcc(c(0.1, 0.2, 0.3, 0.2, 0.1), tau_max = 2)

  # By setting a massive threshold, the only row becomes NA
  res <- pick_optima(mock_wcc, search_method = "global", threshold = 0.99)

  # The summary should print completeness stats but cleanly skip
  # the value distributions since there are no valid lags to summarize
  expect_snapshot(summary(res))
})

test_that("pick_optima respects explicit find_min overrides", {
  # WCC normally looks for peaks (FALSE). Let's force it to look for valleys (TRUE)
  mock_wcc <- create_mock_wcc(c(0.9, 0.5, 0.1, 0.5, 0.9), tau_max = 2)
  res_wcc_min <- pick_optima(
    mock_wcc,
    search_method = "global",
    find_min = TRUE
  )

  expect_equal(res_wcc_min$optimum_value, 0.1)
  expect_equal(res_wcc_min$optimum_lag, 0)

  # WDTW normally looks for valleys (TRUE). Let's force it to look for peaks (FALSE)
  mock_wdtw <- create_mock_wdtw(c(1.1, 5.0, 9.9, 5.0, 1.1), tau_max = 2)
  res_wdtw_max <- pick_optima(
    mock_wdtw,
    search_method = "global",
    find_min = FALSE
  )

  expect_equal(res_wdtw_max$optimum_value, 9.9)
  expect_equal(res_wdtw_max$optimum_lag, 0)
})

test_that("print_optima handles exactly 1 remaining row for pluralization logic", {
  # Exactly 6 rows with default print n=5 leaves exactly 1 remaining row
  df_six <- data.frame(i = 1:6, tau = rep(0, 6), wcc = rep(0.9, 6))
  mock_six <- list(results_df = df_six, settings = list(lag_max = 2))
  class(mock_six) <- c("wcc_res", "list")

  res_six <- pick_optima(mock_six, search_method = "global")

  # This snapshot should capture the singular "... with 1 more row" string
  expect_snapshot(print(res_six))
})

test_that("pick_optima correctly assigns default search methods when NULL", {
  # WCC defaults to 'local' but still requires L_size
  mock_wcc <- create_mock_wcc(c(0.1, 0.5, 0.9, 0.5, 0.1), tau_max = 2)
  res_wcc_default <- pick_optima(mock_wcc, L_size = 3)

  expect_equal(attr(res_wcc_default, "search_method"), "local")
  expect_equal(res_wcc_default$optimum_lag, 0)

  # WDTW defaults to 'global'
  mock_wdtw <- create_mock_wdtw(c(10.0, 5.0, 1.1, 5.0, 10.0), tau_max = 2)
  res_wdtw_default <- pick_optima(mock_wdtw)

  expect_equal(attr(res_wdtw_default, "search_method"), "global")
  expect_equal(res_wdtw_default$optimum_lag, 0)
})

test_that("pick_optima_cpp catches length mismatches", {
  df_corrupted <- data.frame(
    i = rep(100, 4),
    tau = -2:1,
    wcc = c(0.1, 0.5, 0.9, 0.5)
  )

  mock_corrupted <- list(
    results_df = df_corrupted,
    settings = list(lag_max = 2)
  )
  class(mock_corrupted) <- c("wcc_res", "list")

  expect_error(
    pick_optima(mock_corrupted, L_size = 3, search_method = "local"),
    "does not match the expected length"
  )
})

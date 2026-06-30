# Tests for M6: synchrony_multiverse() and suggest_wcc_params() ---------------

# Shared tiny test series (fast: 200 samples at 10 Hz)
make_test_series <- function(n = 200, seed = 42) {
  set.seed(seed)
  t <- seq(0, n / 10, length.out = n)
  list(
    x = sin(2 * pi * 0.5 * t) + rnorm(n, sd = 0.1),
    y = sin(2 * pi * 0.5 * (t - 0.2)) + rnorm(n, sd = 0.1)
  )
}


# =============================================================================
# AC1: synchrony_multiverse() structure + seconds-to-samples conversion -------
# =============================================================================

test_that("AC1: synchrony_multiverse returns bsync_multiverse with correct structure", {
  s <- make_test_series()
  set.seed(1)
  res <- synchrony_multiverse(
    s$x, s$y,
    estimator = "wcc",
    sample_rate = 10,
    window_sec = c(2, 3),
    lag_sec = 1,
    n_surrogates = 19L
  )

  expect_s3_class(res, "bsync_multiverse")
  expect_true(is.list(res))
  expect_named(res, c("grid", "settings", "robustness"))

  # grid is a tibble with expected columns
  expect_s3_class(res$grid, "tbl_df")
  expect_true(all(c(
    "window_sec", "lag_sec", "increment_pct",
    "surrogate_method", "statistic",
    "window_size", "lag_max", "window_increment",
    "n_windows", "observed", "null_mean", "null_sd",
    "es", "p"
  ) %in% names(res$grid)))

  # seconds → samples: window_sec = 2 at 10 Hz → window_size = 20
  row2 <- res$grid[res$grid$window_sec == 2, ][1, ]
  expect_equal(row2$window_size, 20L)

  # lag_max capped at floor(window_size / 2)
  expect_true(all(res$grid$lag_max <= floor(res$grid$window_size / 2L), na.rm = TRUE))
})

test_that("AC1: Invariant 7 — bsync_multiverse carries no raw input or surrogate draws", {
  s <- make_test_series()
  set.seed(2)
  res <- synchrony_multiverse(
    s$x, s$y,
    estimator = "wcc", sample_rate = 10,
    window_sec = 2, lag_sec = 0.5, n_surrogates = 9L
  )
  # No raw data stored
  nms <- names(res)
  expect_false(any(nms %in% c("x", "y", "data")))
  # No raw surrogate matrices
  has_surr_mat <- vapply(res, function(v) is.matrix(v) && nrow(v) == 200, logical(1))
  expect_false(any(has_surr_mat))
})

test_that("AC1: lag_max cap is applied and reported per cell", {
  s <- make_test_series()
  set.seed(3)
  # lag_sec = 5 > window_sec/2 = 1.5 → should be capped
  res <- synchrony_multiverse(
    s$x, s$y,
    estimator = "wcc", sample_rate = 10,
    window_sec = 3, lag_sec = 5, n_surrogates = 9L
  )
  # lag_max should be at most floor(30/2) = 15
  expect_true(all(res$grid$lag_max <= 15L, na.rm = TRUE))
})


# =============================================================================
# AC2: Invariant 2 — matched-null and tail/ES polarity per estimator ----------
# =============================================================================

# Cross-path Invariant 2: when the sole surrogate IS y itself, the aggregate
# computed on the surrogate must equal the aggregate computed on the observed
# data (same statistic), so null_mean == observed. This exercises each
# estimator's multiverse cell helper directly, proving the surrogate path uses
# the matching statistic.

test_that("AC2: WCC — matched null (sole surrogate = y) => null_mean == observed", {
  s <- make_test_series(n = 200)
  y_surr <- matrix(s$y, ncol = 1)
  for (stat in c("mean_abs_z", "peak")) {
    cell <- bsync:::.mv_wcc_cell(
      s$x, s$y,
      window_size = 20L, lag_max = 5L, window_increment = 4L,
      statistic = stat, na.rm = TRUE, y_surr = y_surr
    )
    expect_equal(cell$null_mean, cell$observed,
      tolerance = 1e-9, info = paste("statistic:", stat)
    )
  }
})

test_that("AC2: WDTW — matched null (sole surrogate = y) => null_mean == observed", {
  s <- make_test_series(n = 200)
  y_surr <- matrix(s$y, ncol = 1)
  cell <- bsync:::.mv_wdtw_cell(
    s$x, s$y,
    window_size = 20L, lag_max = 5L, window_increment = 4L,
    scale_method = "global", distance_metric = "L2", y_surr = y_surr
  )
  expect_equal(cell$null_mean, cell$observed, tolerance = 1e-9)
})

test_that("AC2: Granger — matched null (sole surrogate = y) => null_mean == observed both ways", {
  s <- make_test_series(n = 200)
  y_surr <- matrix(s$y, ncol = 1)
  cell <- bsync:::.mv_granger_cell(
    s$x, s$y,
    window_size = 20L, window_increment = 4L, ar_order = 1L, y_surr = y_surr
  )
  expect_equal(cell$null_mean, cell$observed, tolerance = 1e-9) # x -> y
  expect_equal(cell$null_mean_yx, cell$observed_yx, tolerance = 1e-9) # y -> x
})

test_that("AC2: WDTW — ES uses lower-tail polarity (null_mean - obs) / null_sd", {
  s <- make_test_series()
  set.seed(5)
  res <- synchrony_multiverse(
    s$x, s$y,
    estimator = "wdtw", sample_rate = 10,
    window_sec = 2, lag_sec = 0.5, n_surrogates = 19L
  )
  row1 <- res$grid[1, ]
  # Verify ES polarity: for WDTW, ES > 0 means observed < null_mean (better alignment)
  expected_es <- (row1$null_mean - row1$observed) / row1$null_sd
  expect_equal(row1$es, expected_es, tolerance = 1e-10)
})

test_that("AC2: Granger — ES uses upper-tail polarity for x→y", {
  s <- make_test_series()
  set.seed(6)
  res <- synchrony_multiverse(
    s$x, s$y,
    estimator = "wgranger", sample_rate = 10,
    window_sec = 2, n_surrogates = 19L
  )
  row1 <- res$grid[1, ]
  expected_es <- (row1$observed - row1$null_mean) / row1$null_sd
  expect_equal(row1$es, expected_es, tolerance = 1e-10)
})


# =============================================================================
# AC3: Surrogate-reuse efficiency seam ----------------------------------------
# =============================================================================

test_that("AC3: generate_surrogate_phase called once per method, not per cell", {
  s <- make_test_series()
  call_count <- 0L
  # Capture real implementation before mocking to avoid infinite recursion
  real_phase <- bsync:::generate_surrogate_phase
  local_mocked_bindings(
    generate_surrogate_phase = function(y, n_surrogates, trim_odd = FALSE) {
      call_count <<- call_count + 1L
      real_phase(y, n_surrogates = n_surrogates, trim_odd = trim_odd)
    },
    .package = "bsync"
  )
  set.seed(7)
  synchrony_multiverse(
    s$x, s$y,
    estimator = "wcc",
    sample_rate = 10,
    window_sec = c(2, 3, 4), # 3 cells, all sharing "phase"
    lag_sec = 0.5,
    surrogate_method = "phase",
    n_surrogates = 9L
  )
  # Should be called exactly once for "phase", not 3 times (one per cell)
  expect_equal(call_count, 1L)
})

test_that("AC3: Two methods generate two matrices, not 2 × n_cells", {
  s <- make_test_series()
  phase_calls <- 0L
  circ_calls <- 0L
  real_phase <- bsync:::generate_surrogate_phase
  real_circ <- bsync:::generate_surrogate_circular
  local_mocked_bindings(
    generate_surrogate_phase = function(y, n_surrogates, trim_odd = FALSE) {
      phase_calls <<- phase_calls + 1L
      real_phase(y, n_surrogates = n_surrogates, trim_odd = trim_odd)
    },
    generate_surrogate_circular = function(y, n_surrogates, lag_max = NULL) {
      circ_calls <<- circ_calls + 1L
      real_circ(y, n_surrogates = n_surrogates, lag_max = lag_max)
    },
    .package = "bsync"
  )
  set.seed(8)
  # 2 window sizes × 2 surrogate methods = 4 cells
  synchrony_multiverse(
    s$x, s$y,
    estimator = "wcc",
    sample_rate = 10,
    window_sec = c(2, 3),
    lag_sec = 0.5,
    surrogate_method = c("phase", "circular"),
    n_surrogates = 9L
  )
  expect_equal(phase_calls, 1L)
  expect_equal(circ_calls, 1L)
})


# =============================================================================
# AC4: WDTW and Granger adapters produce correct grids ------------------------
# =============================================================================

test_that("AC4: Granger grid has no lag axis; has es_xy/p_xy/es_yx/p_yx columns", {
  s <- make_test_series()
  set.seed(9)
  res <- synchrony_multiverse(
    s$x, s$y,
    estimator = "wgranger", sample_rate = 10,
    window_sec = c(2, 3), n_surrogates = 9L
  )
  expect_true(all(is.na(res$grid$lag_sec)))
  expect_true(all(is.na(res$grid$lag_max)))
  expect_true(all(c("es_yx", "p_yx") %in% names(res$grid)))
  expect_equal(nrow(res$grid), 2L) # 2 window_sec values
})

test_that("AC4: WDTW smoke test — grid has correct columns, n_windows > 0", {
  s <- make_test_series()
  set.seed(10)
  res <- synchrony_multiverse(
    s$x, s$y,
    estimator = "wdtw", sample_rate = 10,
    window_sec = 2, lag_sec = 0.5, n_surrogates = 9L
  )
  expect_s3_class(res, "bsync_multiverse")
  expect_true(res$grid$n_windows[1] > 0L)
  expect_type(res$grid$es[1], "double")
})


# =============================================================================
# AC5: S3 methods for bsync_multiverse ----------------------------------------
# =============================================================================

test_that("AC5: print, summary, tidy, glance, as_tibble work on bsync_multiverse", {
  s <- make_test_series()
  set.seed(11)
  res <- synchrony_multiverse(
    s$x, s$y,
    estimator = "wcc", sample_rate = 10,
    window_sec = 2, lag_sec = 0.5, n_surrogates = 9L
  )

  expect_invisible(print(res))
  expect_invisible(summary(res))

  td <- tidy(res)
  expect_s3_class(td, "tbl_df")
  expect_equal(nrow(td), nrow(res$grid))

  gl <- glance(res)
  expect_s3_class(gl, "tbl_df")
  expect_equal(nrow(gl), 1L)
  expect_true("pct_significant" %in% names(gl))

  tb <- tibble::as_tibble(res)
  expect_equal(tb, td)
})

test_that("AC5: plot.bsync_multiverse runs without error (vdiffr snapshot)", {
  s <- make_test_series()
  set.seed(12)
  res <- synchrony_multiverse(
    s$x, s$y,
    estimator = "wcc", sample_rate = 10,
    window_sec = c(2, 3, 4),
    lag_sec = c(0.5, 1),
    statistic = c("mean_abs_z", "peak"),
    n_surrogates = 9L
  )
  vdiffr::expect_doppelganger(
    "multiverse-spec-curve",
    plot(res)
  )
})


# =============================================================================
# AC6: suggest_wcc_params() reworked signature --------------------------------
# =============================================================================

test_that("AC6: suggest_wcc_params takes x, y, sample_rate and runs PSD", {
  s <- make_test_series(n = 400)
  params <- suggest_wcc_params(
    x           = s$x,
    y           = s$y,
    sample_rate = 10
  )
  expect_named(params, c("window_size", "lag_max", "window_increment", "lag_increment"))
  expect_true(params$window_size > 0)
  expect_true(params$lag_max > 0)
  expect_true(params$lag_max <= floor(params$window_size / 2))
})

test_that("AC6: suggest_wcc_params event_duration_sec override bypasses PSD", {
  s <- make_test_series(n = 400)
  params_psd <- suggest_wcc_params(s$x, s$y, sample_rate = 10)
  params_override <- suggest_wcc_params(
    s$x, s$y,
    sample_rate = 10, event_duration_sec = 2
  )
  # Override of 2s should give window = round(2 * 4 * 10) = 80
  expect_equal(params_override$window_size, 80L)
  # The two paths can differ
  expect_false(identical(params_psd, params_override))
})

test_that("AC6: suggest_wcc_params enforces hard constraints", {
  # Short series to trigger the series-length ceiling
  set.seed(1)
  x <- rnorm(60)
  y <- rnorm(60)
  # event_duration_sec = 10s @ 10Hz → window = 400 > 60/2 = 30, must be capped
  # Both window cap and lag cap warnings fire; suppress and check resulting values
  params <- suppressWarnings(suppressMessages(suggest_wcc_params(
    x, y,
    sample_rate = 10, event_duration_sec = 10
  )))
  expect_true(params$window_size <= 30L)
  expect_true(params$lag_max <= floor(params$window_size / 2))
})

test_that("AC6: suggest_wcc_params errors on bad inputs", {
  set.seed(1)
  x <- rnorm(50)
  y <- rnorm(50)
  expect_error(suggest_wcc_params(x, y, sample_rate = -1), "positive")
  expect_error(suggest_wcc_params(x, y, sample_rate = 10, overlap_pct = 1.5), "\\[0, 1\\)")
})


# =============================================================================
# AC1 robustness: skipped cells reported as NA --------------------------------
# =============================================================================

test_that("Cells too short to compute appear as NA in grid", {
  s <- make_test_series(n = 100)
  set.seed(13)
  # window_sec = 20 s @ 10 Hz → 200 samples, but series is only 100 → too short
  res <- suppressWarnings(synchrony_multiverse(
    s$x, s$y,
    estimator = "wcc", sample_rate = 10,
    window_sec = c(1, 20), # 1s = ok, 20s = too short
    lag_sec = 0.5,
    n_surrogates = 9L
  ))
  expect_equal(nrow(res$grid), 2L)
  # The too-short cell should have NA es
  expect_true(is.na(res$grid$es[res$grid$window_sec == 20]))
  expect_false(is.na(res$grid$es[res$grid$window_sec == 1]))
})

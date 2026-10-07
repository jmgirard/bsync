# Tests for the Granger direction of bsync_multiverse results ------------------
#
# The fixture has x driving y (y depends on x two samples back) and no
# y -> x path, so the two directions give different significance counts,
# median ES, cell rankings, and at least one cell that is significant in one
# direction only (checked in the first test, so the other tests can rely on
# it).

granger_fixture_env <- new.env(parent = emptyenv())

granger_fixture <- function() {
  if (is.null(granger_fixture_env$res)) {
    set.seed(42)
    n <- 300
    x <- as.numeric(stats::filter(rnorm(n), 0.5, method = "recursive"))
    y <- numeric(n)
    for (i in 3:n) {
      y[i] <- 0.3 * y[i - 1] + 0.6 * x[i - 2] + rnorm(1, sd = 0.5)
    }
    set.seed(1)
    granger_fixture_env$res <- synchrony_multiverse(
      x,
      y,
      estimator = "wgranger",
      sample_rate = 10,
      window_sec = c(1, 2, 3, 4),
      increment_pct = c(0.1, 0.25),
      n_surrogates = 39L
    )
  }
  granger_fixture_env$res
}

small_multiverse <- function(estimator) {
  set.seed(3)
  n <- 200
  t <- seq(0, n / 10, length.out = n)
  x <- sin(2 * pi * 0.5 * t) + rnorm(n, sd = 0.1)
  y <- sin(2 * pi * 0.5 * (t - 0.2)) + rnorm(n, sd = 0.1)
  synchrony_multiverse(
    x,
    y,
    estimator = estimator,
    sample_rate = 10,
    window_sec = 2,
    lag_sec = 0.5,
    n_surrogates = 20L
  )
}

# Robustness fields computed here from grid columns, independently of
# multiverse_robustness(). A cell is valid when its ES is not NA, and
# significant when p <= .05.
expected_robustness <- function(es, p) {
  valid <- !is.na(es)
  v_es <- es[valid]
  v_sig <- !is.na(p[valid]) & p[valid] <= 0.05
  list(
    n_cells = length(es),
    n_valid = length(v_es),
    n_significant = sum(v_sig),
    pct_significant = sum(v_sig) / length(v_es),
    median_es = stats::median(v_es),
    iqr_es = stats::IQR(v_es),
    sign_consistent = if (any(v_sig)) mean(v_es[v_sig] > 0) else NA_real_
  )
}

robustness_fields <- c(
  "n_cells",
  "n_valid",
  "n_significant",
  "pct_significant",
  "median_es",
  "iqr_es",
  "sign_consistent"
)


# $robustness_yx ---------------------------------------------------------------

test_that("the fixture separates the two directions", {
  g <- granger_fixture()$grid
  sig_xy <- g$p <= 0.05
  sig_yx <- g$p_yx <= 0.05
  expect_false(sum(sig_xy) == sum(sig_yx))
  expect_false(stats::median(g$es) == stats::median(g$es_yx))
  expect_false(identical(order(g$es), order(g$es_yx)))
  expect_true(any(sig_xy != sig_yx))
})

test_that("a Granger result carries $robustness_yx built from es_yx and p_yx", {
  res <- granger_fixture()
  expect_named(res$robustness, robustness_fields)
  expect_named(res$robustness_yx, robustness_fields)

  expect_equal(
    res$robustness_yx,
    expected_robustness(res$grid$es_yx, res$grid$p_yx)
  )
  expect_equal(res$robustness, expected_robustness(res$grid$es, res$grid$p))
})

test_that("WCC and WDTW results have no $robustness_yx", {
  expect_false("robustness_yx" %in% names(small_multiverse("wcc")))
  expect_false("robustness_yx" %in% names(small_multiverse("wdtw")))
})

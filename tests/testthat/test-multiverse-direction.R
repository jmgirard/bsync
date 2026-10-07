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

expect_yx_refused <- function(expr, estimator) {
  expect_error(
    expr,
    paste0(
      "needs a Granger result[\\s\\S]*estimator is \"",
      estimator,
      "\""
    ),
    class = "rlang_error",
    perl = TRUE
  )
}

expect_bad_direction <- function(expr) {
  expect_error(
    expr,
    "`direction` must be one of \"xy\" or \"yx\", not \"zz\"",
    class = "rlang_error"
  )
}

# print() and summary() write cli messages, one per line. This checks that one
# of the lines contains `text` and drops the others.
expect_line <- function(f, text) {
  suppressMessages(expect_message(f(), text, fixed = TRUE))
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


# glance() ---------------------------------------------------------------------

summary_cols <- c(
  "n_cells",
  "n_valid",
  "n_significant",
  "pct_significant",
  "median_es",
  "iqr_es",
  "sign_consistent"
)

test_that("glance() reports the chosen Granger direction", {
  res <- granger_fixture()

  gl_xy <- glance(res)
  expect_identical(gl_xy$direction, "xy")
  expect_equal(as.list(gl_xy[summary_cols]), res$robustness)

  gl_yx <- glance(res, direction = "yx")
  expect_identical(gl_yx$direction, "yx")
  expect_equal(as.list(gl_yx[summary_cols]), res$robustness_yx)
  expect_false(gl_xy$n_significant == gl_yx$n_significant)
})

test_that("glance() columns of WCC and WDTW results do not change", {
  cols <- c(
    "estimator",
    "n_cells",
    "n_valid",
    "n_significant",
    "pct_significant",
    "median_es",
    "iqr_es",
    "sign_consistent",
    "n_surrogates"
  )
  expect_named(glance(small_multiverse("wcc")), cols)
  expect_named(glance(small_multiverse("wdtw")), cols)
})

test_that("glance() refuses direction = 'yx' for one-direction estimators", {
  expect_yx_refused(glance(small_multiverse("wcc"), direction = "yx"), "wcc")
  expect_yx_refused(glance(small_multiverse("wdtw"), direction = "yx"), "wdtw")
  expect_bad_direction(glance(granger_fixture(), direction = "zz"))
})


# print() and summary() --------------------------------------------------------

sig_line <- function(rb) {
  paste0(
    "Significant (p <= .05): ",
    rb$n_significant,
    " of ",
    rb$n_valid
  )
}

range_line <- function(es) {
  r <- round(range(es, na.rm = TRUE), 3)
  paste0("ES range: [", r[1], ", ", r[2], "]")
}

test_that("print() names the Granger direction and reports its counts", {
  res <- granger_fixture()
  expect_line(function() print(res), "Direction: x -> y")
  expect_line(function() print(res), sig_line(res$robustness))
  expect_line(function() print(res, direction = "yx"), "Direction: y -> x")
  expect_line(
    function() print(res, direction = "yx"),
    sig_line(res$robustness_yx)
  )
  expect_invisible(suppressMessages(print(res, direction = "yx")))
})

test_that("summary() reports the es_yx range and NA count for y -> x", {
  res <- granger_fixture()
  expect_line(function() summary(res), "Direction: x -> y")
  expect_line(function() summary(res), range_line(res$grid$es))

  f <- function() summary(res, direction = "yx")
  expect_line(f, "Direction: y -> x")
  expect_line(f, sig_line(res$robustness_yx))
  expect_line(f, range_line(res$grid$es_yx))
  expect_false(range_line(res$grid$es_yx) == range_line(res$grid$es))

  # The skipped/NA count reads the chosen direction's ES column.
  res_na <- res
  res_na$grid$es_yx[1:3] <- NA_real_
  expect_line(
    function() summary(res_na, direction = "yx"),
    "(including 3 skipped/NA)"
  )
  expect_line(function() summary(res_na), "(including 0 skipped/NA)")
})

test_that("print() and summary() of a WCC result name no direction", {
  res <- small_multiverse("wcc")
  msgs <- testthat::capture_messages(summary(res))
  expect_false(any(grepl("Direction", msgs, fixed = TRUE)))
})

test_that("print() and summary() refuse direction = 'yx' for WCC", {
  res <- small_multiverse("wcc")
  expect_yx_refused(print(res, direction = "yx"), "wcc")
  expect_yx_refused(summary(res, direction = "yx"), "wcc")
  expect_bad_direction(print(granger_fixture(), direction = "zz"))
  expect_bad_direction(summary(granger_fixture(), direction = "zz"))
})


# plot() -----------------------------------------------------------------------

test_that("plot data rank cells by es_yx and mark significance from p_yx", {
  res <- granger_fixture()
  g <- res$grid

  gd <- multiverse_plot_data(res, "yx")
  expect_equal(gd$window_sec, g$window_sec[order(g$es_yx)])
  expect_equal(gd$increment_pct, g$increment_pct[order(g$es_yx)])
  expect_equal(gd$es, sort(g$es_yx))
  expect_equal(
    gd$significant,
    !is.na(g$p_yx[order(g$es_yx)]) & g$p_yx[order(g$es_yx)] <= 0.05
  )

  gd_xy <- multiverse_plot_data(res)
  expect_equal(gd_xy$es, sort(g$es))
  expect_equal(
    gd_xy$significant,
    !is.na(g$p[order(g$es)]) & g$p[order(g$es)] <= 0.05
  )
})

test_that("the plot title names the Granger direction only", {
  res <- granger_fixture()
  rb <- res$robustness
  rb_yx <- res$robustness_yx
  expect_identical(
    multiverse_plot_title(rb, "x -> y"),
    paste0(
      "Synchrony Multiverse (x -> y)  --  ",
      rb$n_significant,
      " / ",
      rb$n_valid,
      " cells significant  |  Median ES = ",
      round(rb$median_es, 2)
    )
  )
  expect_identical(
    multiverse_plot_title(rb_yx, "y -> x"),
    paste0(
      "Synchrony Multiverse (y -> x)  --  ",
      rb_yx$n_significant,
      " / ",
      rb_yx$n_valid,
      " cells significant  |  Median ES = ",
      round(rb_yx$median_es, 2)
    )
  )

  # The plot passes the chosen direction's robustness and label.
  expect_identical(multiverse_direction(res, "yx")$robustness, rb_yx)
  expect_identical(multiverse_direction(res, "yx")$label, "y -> x")
  expect_identical(multiverse_direction(res)$label, "x -> y")

  # A WCC result keeps the title it had before the direction argument.
  wcc <- small_multiverse("wcc")
  wcc_dir <- multiverse_direction(wcc)
  expect_null(wcc_dir$label)
  expect_identical(
    multiverse_plot_title(wcc_dir$robustness, wcc_dir$label),
    paste0(
      "Synchrony Multiverse  --  ",
      wcc$robustness$n_significant,
      " / ",
      wcc$robustness$n_valid,
      " cells significant  |  Median ES = ",
      round(wcc$robustness$median_es, 2)
    )
  )
})

test_that("plot() draws the y -> x specification curve (vdiffr snapshot)", {
  res <- granger_fixture()
  vdiffr::expect_doppelganger(
    "multiverse-granger-yx",
    function() plot(res, direction = "yx")
  )
})

test_that("plot() refuses direction = 'yx' for WCC and unknown directions", {
  expect_yx_refused(plot(small_multiverse("wcc"), direction = "yx"), "wcc")
  expect_bad_direction(plot(granger_fixture(), direction = "zz"))
})

test_that("plot() stops when the chosen direction has no valid cells", {
  res <- granger_fixture()

  # es is fine, es_yx is all NA: only y -> x has nothing to plot.
  no_yx <- res
  no_yx$grid$es_yx <- NA_real_
  expect_error(
    plot(no_yx, direction = "yx"),
    "No valid cells to plot \\(all ES values are NA\\)",
    class = "rlang_error"
  )
  expect_s3_class(multiverse_plot_data(no_yx), "tbl_df")

  # es_yx is fine, es is all NA: only x -> y has nothing to plot.
  no_xy <- res
  no_xy$grid$es <- NA_real_
  expect_error(
    plot(no_xy, direction = "xy"),
    "No valid cells to plot \\(all ES values are NA\\)",
    class = "rlang_error"
  )
  expect_s3_class(multiverse_plot_data(no_xy, "yx"), "tbl_df")
})


# Granger results without $robustness_yx or the _yx columns ---------------------

test_that("y -> x is computed from the grid when $robustness_yx is missing", {
  # A Granger result saved before $robustness_yx existed. Some p_yx values
  # are set significant, so sign_consistent is a number, not NA.
  old <- granger_fixture()
  old$robustness_yx <- NULL
  old$grid$p_yx[c(1, 4, 6)] <- 0.01
  expected <- expected_robustness(old$grid$es_yx, old$grid$p_yx)
  expect_equal(expected$n_significant, 3L)
  expect_false(is.na(expected$sign_consistent))

  expect_equal(multiverse_direction(old, "yx")$robustness, expected)
  gl <- glance(old, direction = "yx")
  expect_named(gl, c("estimator", "direction", summary_cols, "n_surrogates"))
  expect_equal(as.list(gl[summary_cols]), expected)
  expect_line(function() print(old, direction = "yx"), sig_line(expected))
  expect_identical(
    multiverse_plot_title(multiverse_direction(old, "yx")$robustness),
    multiverse_plot_title(expected)
  )
})

test_that("y -> x stops with a message when the grid has no _yx columns", {
  res <- granger_fixture()
  res$grid$es_yx <- NULL
  res$grid$p_yx <- NULL
  msg <- "needs the es_yx and p_yx grid columns"
  expect_error(glance(res, direction = "yx"), msg, class = "rlang_error")
  expect_error(summary(res, direction = "yx"), msg, class = "rlang_error")
  expect_error(plot(res, direction = "yx"), msg, class = "rlang_error")
})

test_that("the refusal names an unknown estimator when settings lack one", {
  res <- small_multiverse("wcc")
  res$settings$estimator <- NULL
  expect_error(
    glance(res, direction = "yx"),
    "estimator is \"unknown\"",
    class = "rlang_error"
  )
})

test_that("summary() says 'none' when the direction has no computable ES", {
  res <- granger_fixture()
  res$grid$es_yx <- NA_real_
  f <- function() summary(res, direction = "yx")
  expect_no_warning(suppressMessages(f()))
  expect_line(f, "ES range: none (no computable ES)")
  expect_invisible(suppressMessages(f()))
})

test_that("summary(direction = 'yx') returns its input invisibly", {
  res <- granger_fixture()
  expect_invisible(suppressMessages(summary(res, direction = "yx")))
})

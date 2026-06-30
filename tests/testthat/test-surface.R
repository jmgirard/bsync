# Tests for the shared windowed-surface infrastructure (M5) ------------------
#
# AC1 characterization: results_df and aggregate values are bit-identical to
# the pre-refactor implementation on sim_dyad (captured 2026-06-29 before any
# M5 code changes — see M5 plan in CLAUDE.md). Any regression shows up here.
#
# AC2 validator: shared validation path aborts with the same messages as before.
#
# AC3 superclass: all three result objects inherit "bsync_surface" and carry
# a named-numeric $aggregate slot.


# AC1: characterization (output-preserving refactor) --------------------------

test_that("AC1: wcc results_df and aggregate are bit-identical after refactor", {
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B

  res <- wcc(x, y, window_size = 96, lag_max = 10)

  expect_equal(nrow(res$results_df), 47985L)
  expect_equal(names(res$results_df), c("i", "tau", "wcc"))
  expect_equal(res$aggregate[[1]], 0.090014146265991, tolerance = 1e-12)

  res_peak <- wcc(x, y, window_size = 96, lag_max = 10, statistic = "peak")
  expect_equal(res_peak$aggregate[[1]], 0.236893791365889, tolerance = 1e-12)
})

test_that("AC1: wdtw results_df and aggregate are bit-identical after refactor", {
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B

  res <- wdtw(x, y, window_size = 96, lag_max = 10)

  expect_equal(nrow(res$results_df), 47985L)
  expect_equal(names(res$results_df), c("i", "tau", "dtw_dist"))
  expect_equal(res$aggregate[["mean_distance"]], 51.400395173507370, tolerance = 1e-10)
})

test_that("AC1: wgranger results_df and aggregate are bit-identical after refactor", {
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B

  res <- wgranger(x, y, window_size = 96, ar_order = 1)

  expect_equal(nrow(res$results_df), 2305L)
  expect_equal(names(res$results_df), c("i", "f_xy", "p_xy", "f_yx", "p_yx"))
  expect_equal(res$aggregate[["f_xy"]], 1.092260756425073, tolerance = 1e-10)
  expect_equal(res$aggregate[["f_yx"]], 1.848471978705106, tolerance = 1e-10)
})


# AC1: one grid builder (grep is the authoritative check; this tests behavior) -

test_that("AC1: build_surface_grid matches create_*_df grid math (lagged)", {
  # Confirm grid builder produces the same (i, tau) pairs as the estimators.
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B

  res_wcc <- wcc(x, y, window_size = 96, lag_max = 10)
  grid <- build_surface_grid(
    n_x = length(x), window_size = 96, window_increment = 1,
    lag_max = 10, lag_increment = 1, lagged = TRUE
  )

  expect_equal(grid$i_vals, res_wcc$results_df$i)
  expect_equal(grid$tau_vals, res_wcc$results_df$tau)
  expect_equal(grid$w_max, 95L)
  expect_equal(grid$n_r, 2285L)
})

test_that("AC1: build_surface_grid matches create_wgranger_df grid math (lag-free)", {
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B

  res_wg <- wgranger(x, y, window_size = 96)
  grid <- build_surface_grid(
    n_x = length(x), window_size = 96, window_increment = 1,
    lagged = FALSE
  )

  expect_equal(grid$i_vals, res_wg$results_df$i)
  expect_equal(grid$w_max, 95L)
  expect_equal(grid$n_r, 2305L)
})

test_that("AC1: build_surface_grid aborts for too-short series", {
  expect_error(
    build_surface_grid(n_x = 5, window_size = 10, lag_max = 3, lagged = TRUE),
    "too short"
  )
  expect_error(
    build_surface_grid(n_x = 5, window_size = 10, lagged = FALSE),
    "too short"
  )
})


# AC2: shared validator -------------------------------------------------------

test_that("AC2: validate_series catches type and length errors", {
  expect_error(validate_series(c("a"), 1:5), "numeric vector")
  expect_error(validate_series(1:5, c("a")), "numeric vector")
  expect_error(validate_series(1:3, 1:5), "same length")
  expect_error(validate_series(1:5, 1:5, time = c("t")), "numeric vector")
  expect_error(validate_series(1:5, 1:5, time = 1:3), "same length")
})

test_that("AC2: validate_window_params catches bad params", {
  expect_error(validate_window_params(window_size = -1), "positive integer")
  expect_error(validate_window_params(window_size = 5, window_increment = 0), "positive integer")
  expect_error(validate_window_params(window_size = 5, lag_max = 0), "positive integer")
  expect_error(validate_window_params(window_size = 5, lag_increment = -1), "positive integer")
  expect_error(validate_window_params(window_size = 5, ar_order = 0), "positive integer")
  expect_error(validate_window_params(window_size = 5, na.rm = "yes"), "logical value")
})


# AC3: bsync_surface superclass + $aggregate slot -----------------------------

test_that("AC3: all three result objects inherit bsync_surface", {
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B

  wcc_res <- wcc(x, y, window_size = 96, lag_max = 10)
  wdtw_res <- wdtw(x, y, window_size = 96, lag_max = 10)
  wg_res <- wgranger(x, y, window_size = 96)

  expect_true(inherits(wcc_res, "bsync_surface"))
  expect_true(inherits(wdtw_res, "bsync_surface"))
  expect_true(inherits(wg_res, "bsync_surface"))

  # Leaf classes still present
  expect_s3_class(wcc_res, "wcc_res")
  expect_s3_class(wdtw_res, "wdtw_res")
  expect_s3_class(wg_res, "wgranger_res")
})

test_that("AC3: $aggregate is a named numeric on all three", {
  x <- sim_dyad$x_A
  y <- sim_dyad$x_B

  wcc_res <- wcc(x, y, window_size = 96, lag_max = 10)
  wcc_pk <- wcc(x, y, window_size = 96, lag_max = 10, statistic = "peak")
  wdtw_res <- wdtw(x, y, window_size = 96, lag_max = 10)
  wg_res <- wgranger(x, y, window_size = 96)

  expect_named(wcc_res$aggregate, "mean_abs_z")
  expect_named(wcc_pk$aggregate, "peak")
  expect_named(wdtw_res$aggregate, "mean_distance")
  expect_named(wg_res$aggregate, c("f_xy", "f_yx"))
})

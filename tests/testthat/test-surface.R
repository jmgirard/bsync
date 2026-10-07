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


# Shared read-only fixtures ---------------------------------------------------
# A single full-series WDTW surface on sim_dyad costs ~19s (DTW is O(window^2)
# per cell). The AC1 characterization values are frozen on the full series at
# window_size = 96 / lag_max = 10, and the AC3/AC6 structure tests assert
# against that same surface (e.g. n_windows == 2285). So each surface is built
# exactly once here and reused below — identical objects, identical assertions,
# no loss of rigor, but the expensive WDTW surface is computed a single time.
x_ref <- sim_dyad$x_A
y_ref <- sim_dyad$x_B
wcc_ref <- wcc(x_ref, y_ref, window_size = 96, lag_max = 10)
wcc_peak_ref <- wcc(x_ref, y_ref, window_size = 96, lag_max = 10, statistic = "peak")
wdtw_ref <- wdtw(x_ref, y_ref, window_size = 96, lag_max = 10)
wg_ref <- wgranger(x_ref, y_ref, window_size = 96, ar_order = 1)


# AC1: characterization (output-preserving refactor) --------------------------

test_that("AC1: wcc results_df and aggregate are bit-identical after refactor", {
  expect_equal(nrow(wcc_ref$results_df), 47985L)
  expect_equal(names(wcc_ref$results_df), c("i", "tau", "wcc"))
  expect_equal(wcc_ref$aggregate[[1]], 0.090014146265991, tolerance = 1e-12)

  expect_equal(wcc_peak_ref$aggregate[[1]], 0.236893791365889, tolerance = 1e-12)
})

test_that("AC1: wdtw results_df and aggregate are bit-identical after refactor", {
  expect_equal(nrow(wdtw_ref$results_df), 47985L)
  expect_equal(names(wdtw_ref$results_df), c("i", "tau", "dtw_dist"))
  expect_equal(wdtw_ref$aggregate[["mean_distance"]], 51.400395173507370, tolerance = 1e-10)
})

test_that("AC1: wgranger results_df and aggregate are bit-identical after refactor", {
  expect_equal(nrow(wg_ref$results_df), 2305L)
  expect_equal(names(wg_ref$results_df), c("i", "f_xy", "p_xy", "f_yx", "p_yx"))
  expect_equal(wg_ref$aggregate[["f_xy"]], 1.092260756425073, tolerance = 1e-10)
  expect_equal(wg_ref$aggregate[["f_yx"]], 1.848471978705106, tolerance = 1e-10)
})


# AC1: one grid builder (grep is the authoritative check; this tests behavior) -

test_that("AC1: build_surface_grid matches create_*_df grid math (lagged)", {
  # Confirm grid builder produces the same (i, tau) pairs as the estimators.
  grid <- build_surface_grid(
    n_x = length(x_ref), window_size = 96, window_increment = 1,
    lag_max = 10, lag_increment = 1, lagged = TRUE
  )

  expect_equal(grid$i_vals, wcc_ref$results_df$i)
  expect_equal(grid$tau_vals, wcc_ref$results_df$tau)
  expect_equal(grid$w_max, 95L)
  expect_equal(grid$n_r, 2285L)
})

test_that("AC1: build_surface_grid matches create_wgranger_df grid math (lag-free)", {
  grid <- build_surface_grid(
    n_x = length(x_ref), window_size = 96, window_increment = 1,
    lagged = FALSE
  )

  expect_equal(grid$i_vals, wg_ref$results_df$i)
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


# Window count: every window that fits ----------------------------------------

# Window starts found by an explicit loop: start i is kept when its lagged
# window lies inside 1..n_x, stepping by the increment from 1 + lag_max.
loop_starts <- function(n_x, window_size, window_increment, lag_max) {
  starts <- integer(0)
  i <- 1L + lag_max
  while (i + window_size - 1L + lag_max <= n_x) {
    if (i - lag_max >= 1L) starts <- c(starts, i)
    i <- i + window_increment
  }
  starts
}

# Returns a character vector naming each mismatch between the grid and the
# loop (empty when they agree).
sweep_grid <- function(lagged) {
  bad <- character(0)
  lag_range <- if (lagged) 1:3 else 0L
  for (n_x in 1:40) for (ws in 2:8) for (inc in 1:4) for (lag in lag_range) {
    expected <- loop_starts(n_x, ws, inc, lag)
    min_n <- ws + 2L * lag
    label <- sprintf("n_x=%d ws=%d inc=%d lag=%d", n_x, ws, inc, lag)
    res <- tryCatch(
      build_surface_grid(
        n_x = n_x, window_size = ws, window_increment = inc,
        lag_max = if (lagged) lag else NULL, lagged = lagged
      ),
      error = function(e) e
    )
    if (inherits(res, "error")) {
      msg <- conditionMessage(res)
      if (length(expected) > 0L) {
        bad <- c(bad, paste(label, "aborted but", length(expected), "fit"))
      } else if (!grepl("too short", msg, fixed = TRUE) ||
        !grepl(sprintf("at least %d samples", min_n), msg, fixed = TRUE)) {
        bad <- c(bad, paste(label, "wrong message:", msg))
      }
      next
    }
    if (length(expected) == 0L) {
      bad <- c(bad, paste(label, "returned windows but none fit"))
      next
    }
    starts <- unique(res$i_vals)
    if (!identical(as.integer(starts), expected) ||
      res$n_r != length(expected)) {
      bad <- c(bad, paste(label, "starts differ"))
    }
    if (lagged && !identical(
      as.integer(res$tau_vals), rep(-lag:lag, each = length(expected))
    )) {
      bad <- c(bad, paste(label, "tau_vals differ"))
    }
  }
  bad
}

test_that("build_surface_grid() returns every window start that fits (lagged)", {
  expect_identical(sweep_grid(lagged = TRUE), character(0))
})

test_that("build_surface_grid() returns every window start that fits (lag-free)", {
  expect_identical(sweep_grid(lagged = FALSE), character(0))
})

test_that("build_surface_grid() keeps the pre-fix grid when window_increment = 1", {
  # Pre-fix formula: n_r = floor((n_x - w_max - 2 * lag_max) / increment),
  # starts at 1 + lag_max. With increment 1 it must give the same grid.
  bad <- character(0)
  for (lagged in c(TRUE, FALSE)) {
    lag_range <- if (lagged) 1:3 else 0L
    for (n_x in 1:40) for (ws in 2:8) for (lag in lag_range) {
      old_n_r <- floor((n_x - (ws - 1L) - 2L * lag) / 1L)
      res <- tryCatch(
        build_surface_grid(
          n_x = n_x, window_size = ws, window_increment = 1L,
          lag_max = if (lagged) lag else NULL, lagged = lagged
        ),
        error = function(e) NULL
      )
      label <- sprintf("n_x=%d ws=%d lag=%d lagged=%s", n_x, ws, lag, lagged)
      if (old_n_r < 1L) {
        if (!is.null(res)) bad <- c(bad, paste(label, "no longer aborts"))
        next
      }
      if (is.null(res) || res$n_r != old_n_r ||
        !identical(
          as.integer(unique(res$i_vals)),
          as.integer(1L + lag + (seq_len(old_n_r) - 1L))
        )) {
        bad <- c(bad, paste(label, "differs from pre-fix grid"))
      }
    }
  }
  expect_identical(bad, character(0))
})

test_that("build_surface_grid() returns one window when exactly one fits", {
  grid <- build_surface_grid(
    n_x = 5, window_size = 3, window_increment = 5, lag_max = 1
  )
  expect_equal(grid$n_r, 1L)
  expect_equal(unique(grid$i_vals), 2L)
  expect_equal(grid$tau_vals, -1:1)
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
  expect_true(inherits(wcc_ref, "bsync_surface"))
  expect_true(inherits(wdtw_ref, "bsync_surface"))
  expect_true(inherits(wg_ref, "bsync_surface"))

  # Leaf classes still present
  expect_s3_class(wcc_ref, "wcc_res")
  expect_s3_class(wdtw_ref, "wdtw_res")
  expect_s3_class(wg_ref, "wgranger_res")
})

test_that("AC3: $aggregate is a named numeric on all three", {
  expect_named(wcc_ref$aggregate, "mean_abs_z")
  expect_named(wcc_peak_ref$aggregate, "peak")
  expect_named(wdtw_ref$aggregate, "mean_distance")
  expect_named(wg_ref$aggregate, c("f_xy", "f_yx"))
})


# AC6: tidy interface (tidy / glance / as_tibble) -----------------------------

test_that("AC6: tidy.bsync_surface returns tibble matching results_df", {
  t <- generics::tidy(wcc_ref)

  expect_s3_class(t, "tbl_df")
  expect_equal(nrow(t), nrow(wcc_ref$results_df))
  expect_equal(names(t), c("i", "tau", "wcc"))
  expect_equal(t$wcc, wcc_ref$results_df$wcc)
})

test_that("AC6: tidy / as_tibble return the same tibble", {
  expect_equal(generics::tidy(wdtw_ref), tibble::as_tibble(wdtw_ref))
})

test_that("AC6: glance.bsync_surface — WCC has 1 row with correct fields", {
  g <- generics::glance(wcc_ref)

  expect_s3_class(g, "tbl_df")
  expect_equal(nrow(g), 1L)
  expect_true("mean_abs_z" %in% names(g))
  expect_equal(g$mean_abs_z, wcc_ref$aggregate[["mean_abs_z"]], tolerance = 1e-12)
  expect_equal(g$n_windows, 2285L)
  expect_equal(g$window_size, 96L)
  expect_equal(g$lag_max, 10L)
  expect_equal(g$statistic, "mean_abs_z")
})

test_that("AC6: glance.bsync_surface — WDTW has 1 row with correct fields", {
  g <- generics::glance(wdtw_ref)

  expect_s3_class(g, "tbl_df")
  expect_equal(nrow(g), 1L)
  expect_true("mean_distance" %in% names(g))
  expect_equal(g$mean_distance, wdtw_ref$aggregate[["mean_distance"]], tolerance = 1e-10)
  expect_equal(g$n_windows, 2285L)
  expect_true("scale_method" %in% names(g))
})

test_that("AC6: glance.bsync_surface — Granger has f_xy and f_yx columns", {
  g <- generics::glance(wg_ref)

  expect_s3_class(g, "tbl_df")
  expect_equal(nrow(g), 1L)
  expect_true(all(c("f_xy", "f_yx") %in% names(g)))
  expect_equal(g$f_xy, wg_ref$aggregate[["f_xy"]], tolerance = 1e-10)
  expect_equal(g$f_yx, wg_ref$aggregate[["f_yx"]], tolerance = 1e-10)
  expect_equal(g$n_windows, 2305L)
  expect_equal(g$ar_order, 1L)
})

test_that("AC6: glance() outputs from two runs can be bound into a tibble", {
  combined <- dplyr::bind_rows(generics::glance(wcc_ref), generics::glance(wcc_peak_ref))

  expect_equal(nrow(combined), 2L)
  expect_equal(combined$statistic, c("mean_abs_z", "peak"))
})

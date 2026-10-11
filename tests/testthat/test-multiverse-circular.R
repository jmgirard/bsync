# Circular cells of synchrony_multiverse() (D-005): each cell shifts y by at
# least 2 * lag_max of its own lag (wcc, wdtw) or by at least its window
# (wgranger), uses min(n_surrogates, k) distinct shifts of its k valid shifts,
# and is skipped when it fits its window with k < 19.

# Run `fn` with the cell functions wrapped, and return the result with the
# surrogate matrix that each computed cell received.
capture_cells <- function(fn) {
  seen <- list()
  real_wcc <- bsync:::.mv_wcc_cell
  real_wdtw <- bsync:::.mv_wdtw_cell
  real_granger <- bsync:::.mv_granger_cell
  record <- function(window_size, lag_max, y_surr) {
    seen[[length(seen) + 1L]] <<- list(
      window_size = window_size,
      lag_max = lag_max,
      y_surr = y_surr
    )
  }
  local_mocked_bindings(
    .mv_wcc_cell = function(
      x,
      y,
      window_size,
      lag_max,
      window_increment,
      statistic,
      na.rm,
      y_surr
    ) {
      record(window_size, lag_max, y_surr)
      real_wcc(
        x,
        y,
        window_size,
        lag_max,
        window_increment,
        statistic,
        na.rm,
        y_surr
      )
    },
    .mv_wdtw_cell = function(
      x,
      y,
      window_size,
      lag_max,
      window_increment,
      scale_method,
      distance_metric,
      y_surr
    ) {
      record(window_size, lag_max, y_surr)
      real_wdtw(
        x,
        y,
        window_size,
        lag_max,
        window_increment,
        scale_method,
        distance_metric,
        y_surr
      )
    },
    .mv_granger_cell = function(
      x,
      y,
      window_size,
      window_increment,
      ar_order,
      y_surr
    ) {
      record(window_size, NA_integer_, y_surr)
      real_granger(x, y, window_size, window_increment, ar_order, y_surr)
    },
    .package = "bsync"
  )
  res <- fn()
  list(res = res, cells = seen)
}

# The shift of each column of `y_surr`: column j equals y moved back by s
# samples, y[((i + s - 1) %% n) + 1]. `y` must have distinct values. A column
# that is no rotation of y fails the test.
recover_shifts <- function(y, y_surr) {
  n <- length(y)
  vapply(
    seq_len(ncol(y_surr)),
    function(j) {
      s <- match(y_surr[1, j], y) - 1L
      expect_identical(y_surr[, j], y[((seq_len(n) + s - 1L) %% n) + 1L])
      s
    },
    integer(1)
  )
}

# Cells captured for one window and lag (lag NA for wgranger).
cell_shifts <- function(out, y, window_size, lag_max = NA_integer_) {
  hits <- Filter(
    function(cl) {
      cl$window_size == window_size &&
        identical(
          as.integer(cl$lag_max),
          as.integer(lag_max)
        )
    },
    out$cells
  )
  expect_gte(length(hits), 1L)
  lapply(hits, function(cl) recover_shifts(y, cl$y_surr))
}

test_that("wgranger circular cells shift by at least one window (n = 200 to 204)", {
  for (n in 200:204) {
    set.seed(n)
    x <- rnorm(n)
    y <- rnorm(n)
    out <- NULL
    expect_no_error(
      out <- capture_cells(function() {
        suppressMessages(synchrony_multiverse(
          x,
          y,
          estimator = "wgranger",
          sample_rate = 10,
          window_sec = c(2, 5),
          surrogate_method = "circular",
          n_surrogates = 19L
        ))
      })
    )
    for (w in c(20L, 50L)) {
      for (s in cell_shifts(out, y, w)) {
        expect_true(all(s >= w & s <= n - w), label = paste("n", n, "w", w))
        expect_length(unique(s), 19L)
      }
    }
    expect_identical(out$res$grid$n_surrogates, c(19L, 19L))
  }
})

test_that("wcc and wdtw circular cells shift by twice their own lag", {
  set.seed(401)
  x <- rnorm(401)
  y <- rnorm(401)
  for (est in c("wcc", "wdtw")) {
    # Lags 20 and 100 are capped to 5 for window 10, so those two cells repeat
    # the cell (10, 5) and are dropped. Cell (200, 100) has shifts 200 to 201
    # only, k = 2, so it is skipped.
    warns <- list()
    out <- withCallingHandlers(
      capture_cells(function() {
        suppressMessages(synchrony_multiverse(
          x,
          y,
          estimator = est,
          sample_rate = 1,
          window_sec = c(10, 200),
          lag_sec = c(5, 20, 100),
          surrogate_method = "circular",
          n_surrogates = 50L
        ))
      }),
      bsync_circular_skipped = function(w) {
        warns[[length(warns) + 1L]] <<- w
        invokeRestart("muffleWarning")
      }
    )
    s105 <- cell_shifts(out, y, 10L, 5L)[[1]]
    expect_length(unique(s105), 50L)
    expect_true(all(s105 >= 10L & s105 <= 391L))
    # On main, every cell used shifts 200:201, the range of the largest lag.
    expect_true(any(!s105 %in% 200:201))
    s20020 <- cell_shifts(out, y, 200L, 20L)[[1]]
    expect_length(unique(s20020), 50L)
    expect_true(all(s20020 >= 40L & s20020 <= 361L))

    expect_length(warns, 1L)
    expect_identical(warns[[1]]$window_size, 200L)
    expect_identical(warns[[1]]$lag_max, 100L)
    expect_identical(warns[[1]]$n_shifts, 2L)
    g <- out$res$grid
    skipped <- g[g$window_size == 200L & g$lag_max == 100L, ]
    expect_true(is.na(skipped$es))
    expect_identical(skipped$n_surrogates, 0L)
  }
})

test_that("a circular cell uses min(n_surrogates, k) distinct shifts", {
  # wcc and wdtw: n = 49, window 10, lag 5 gives shifts 10 to 39, k = 30.
  set.seed(49)
  x <- rnorm(49)
  y <- rnorm(49)
  for (est in c("wcc", "wdtw")) {
    out <- capture_cells(function() {
      suppressMessages(synchrony_multiverse(
        x,
        y,
        estimator = est,
        sample_rate = 1,
        window_sec = 10,
        lag_sec = 5,
        surrogate_method = "circular",
        n_surrogates = 50L
      ))
    })
    s <- cell_shifts(out, y, 10L, 5L)[[1]]
    expect_identical(sort(s), 10:39)
    expect_identical(out$res$grid$n_surrogates, 30L)
  }
  # wgranger: n = 69, window 20 gives shifts 20 to 49, k = 30.
  set.seed(69)
  x <- rnorm(69)
  y <- rnorm(69)
  out <- capture_cells(function() {
    suppressMessages(synchrony_multiverse(
      x,
      y,
      estimator = "wgranger",
      sample_rate = 1,
      window_sec = 20,
      surrogate_method = "circular",
      n_surrogates = 50L
    ))
  })
  s <- cell_shifts(out, y, 20L)[[1]]
  expect_identical(sort(s), 20:49)
  expect_identical(out$res$grid$n_surrogates, 30L)
})

# Run one circular cell and return the result with the circular-skip
# warnings it gave.
run_circular <- function(
  n,
  estimator,
  window_sec,
  lag_sec = NULL,
  n_surrogates = 19L,
  seed = 1,
  quiet = TRUE
) {
  set.seed(seed)
  x <- rnorm(n)
  y <- rnorm(n)
  warns <- list()
  run <- function() {
    synchrony_multiverse(
      x,
      y,
      estimator = estimator,
      sample_rate = 1,
      window_sec = window_sec,
      lag_sec = lag_sec,
      surrogate_method = "circular",
      n_surrogates = n_surrogates
    )
  }
  res <- withCallingHandlers(
    if (quiet) suppressMessages(run()) else run(),
    bsync_circular_skipped = function(w) {
      warns[[length(warns) + 1L]] <<- w
      invokeRestart("muffleWarning")
    },
    bsync_few_surrogates = function(w) invokeRestart("muffleWarning")
  )
  list(res = res, warns = warns)
}

test_that("a wcc circular cell with 18 shifts is skipped and 19 computes", {
  # Window 4 and lag 1: shifts 2 to n - 2, so k = n - 3.
  r18 <- run_circular(21, "wcc", window_sec = 4, lag_sec = 1)
  expect_length(r18$warns, 1L)
  expect_identical(r18$warns[[1]]$n_shifts, 18L)
  expect_match(conditionMessage(r18$warns[[1]]), "18 distinct\\s+shifts")
  expect_true(is.na(r18$res$grid$es))
  expect_true(is.na(r18$res$grid$p))
  expect_identical(r18$res$grid$n_surrogates, 0L)

  r19 <- run_circular(22, "wcc", window_sec = 4, lag_sec = 1)
  expect_length(r19$warns, 0L)
  expect_false(is.na(r19$res$grid$p))
  expect_identical(r19$res$grid$n_surrogates, 19L)
})

test_that("a wgranger circular cell with 18 shifts is skipped and 19 computes", {
  # Window 10: shifts 10 to n - 10, so k = n - 19.
  r18 <- run_circular(37, "wgranger", window_sec = 10)
  expect_length(r18$warns, 1L)
  expect_identical(r18$warns[[1]]$n_shifts, 18L)
  expect_true(is.na(r18$res$grid$es))
  expect_identical(r18$res$grid$n_surrogates, 0L)

  r19 <- run_circular(38, "wgranger", window_sec = 10)
  expect_length(r19$warns, 0L)
  expect_false(is.na(r19$res$grid$p))
  expect_identical(r19$res$grid$n_surrogates, 19L)
})

test_that("the skip below 19 shifts holds when fewer surrogates are asked", {
  # k = 18 with 9 surrogates asked: skipped. k = 27 with 9 asked: 9 used.
  r18 <- run_circular(21, "wcc", window_sec = 4, lag_sec = 1, n_surrogates = 9L)
  expect_length(r18$warns, 1L)
  expect_identical(r18$res$grid$n_surrogates, 0L)
  r27 <- run_circular(30, "wcc", window_sec = 4, lag_sec = 1, n_surrogates = 9L)
  expect_length(r27$warns, 0L)
  expect_identical(r27$res$grid$n_surrogates, 9L)
})

test_that("a circular cell too short for its window gives no circular warning", {
  # wcc: window 30 needs 32 samples. wgranger: window 40 needs 40 samples.
  for (case in list(
    list(est = "wcc", w = 30, l = 1),
    list(est = "wgranger", w = 40, l = NULL)
  )) {
    expect_message(
      r <- run_circular(
        21,
        case$est,
        window_sec = case$w,
        lag_sec = case$l,
        quiet = FALSE
      ),
      "series too short"
    )
    expect_length(r$warns, 0L)
    expect_identical(r$res$grid$n_surrogates, 0L)
    expect_true(is.na(r$res$grid$es))
  }
})

test_that("a lag-0 circular cell (1-sample window) shifts by 1 to n - 1", {
  # The lag cap makes lag 0 for a 1-sample window. Its minimum shift is 1
  # sample, so its shifts run from 1 to n - 1. Before the floor, the cell drew
  # shift 0 or n and the whole call stopped.
  set.seed(100)
  x <- rnorm(100)
  y <- rnorm(100)
  for (est in c("wcc", "wdtw")) {
    out <- NULL
    expect_no_error(
      out <- capture_cells(function() {
        suppressWarnings(suppressMessages(synchrony_multiverse(
          x,
          y,
          estimator = est,
          sample_rate = 1,
          window_sec = c(1, 10),
          lag_sec = 2,
          surrogate_method = "circular",
          n_surrogates = 99L
        )))
      })
    )
    g <- out$res$grid
    expect_identical(g$lag_max, c(0L, 2L))
    s0 <- cell_shifts(out, y, 1L, 0L)[[1]]
    expect_identical(sort(s0), 1:99)
    # Lag 2: shifts 4 to 96, k = 93, below the 99 requested.
    expect_identical(g$n_surrogates[g$lag_max == 2L], 93L)
    expect_false(is.na(g$p[g$lag_max == 2L]))
  }
})

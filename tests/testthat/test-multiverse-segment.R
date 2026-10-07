library(testthat)

# Segment surrogates in synchrony_multiverse() and autotune_wcc(). Design:
# D-004 and RR01 Q4 (cairn/reviews/archive/RR01-segment-surrogate-size.md).
# A segment cell uses segments of window + 2 * lag samples (the window for
# Granger) and steps one segment at a time. sample_rate = 1 below, so
# seconds equal samples.

ar1 <- function(n, phi = 0.7) {
  as.numeric(stats::filter(stats::rnorm(n), phi, method = "recursive"))
}

seg_pair <- local({
  set.seed(501)
  list(x = ar1(400), y = ar1(400))
})

count_conditions <- function(expr, class) {
  n <- 0L
  withCallingHandlers(
    expr,
    condition = function(cnd) {
      if (inherits(cnd, class)) n <<- n + 1L
    }
  )
  n
}

test_that("a one-cell wcc segment multiverse matches wcc_surrogate()", {
  w <- 16
  lag <- 2
  s <- w + 2 * lag
  set.seed(601)
  mv <- suppressMessages(synchrony_multiverse(
    seg_pair$x,
    seg_pair$y,
    estimator = "wcc",
    sample_rate = 1,
    window_sec = w,
    lag_sec = lag,
    surrogate_method = "segment",
    n_surrogates = 19
  ))
  set.seed(601)
  surr <- suppressMessages(
    generate_surrogate_segment(seg_pair$y, segment_size = s, n_surrogates = 19)
  )
  ref <- wcc_surrogate(
    seg_pair$x,
    seg_pair$y,
    surr,
    window_size = w,
    lag_max = lag,
    window_increment = s
  )
  expect_identical(nrow(mv$grid), 1L)
  expect_equal(mv$grid$p, ref$p_value)
  # Doubles, the grid column type of every method on main.
  expect_identical(mv$grid$window_increment, as.double(s))
  expect_true(is.na(mv$grid$increment_pct))
})

test_that("a one-cell wdtw segment multiverse matches wdtw_surrogate()", {
  w <- 16
  lag <- 2
  s <- w + 2 * lag
  set.seed(602)
  mv <- suppressMessages(synchrony_multiverse(
    seg_pair$x,
    seg_pair$y,
    estimator = "wdtw",
    sample_rate = 1,
    window_sec = w,
    lag_sec = lag,
    surrogate_method = "segment",
    n_surrogates = 19
  ))
  set.seed(602)
  surr <- suppressMessages(
    generate_surrogate_segment(seg_pair$y, segment_size = s, n_surrogates = 19)
  )
  ref <- wdtw_surrogate(
    seg_pair$x,
    seg_pair$y,
    surr,
    window_size = w,
    lag_max = lag,
    window_increment = s
  )
  expect_equal(mv$grid$p, ref$p_value)
  expect_identical(mv$grid$window_increment, as.double(s))
})

test_that("a one-cell wgranger segment multiverse matches wgranger_surrogate()", {
  w <- 32
  set.seed(603)
  mv <- suppressMessages(synchrony_multiverse(
    seg_pair$x,
    seg_pair$y,
    estimator = "wgranger",
    sample_rate = 1,
    window_sec = w,
    surrogate_method = "segment",
    n_surrogates = 19,
    ar_order = 2L
  ))
  set.seed(603)
  surr <- suppressMessages(
    generate_surrogate_segment(seg_pair$y, segment_size = w, n_surrogates = 19)
  )
  ref <- wgranger_surrogate(
    seg_pair$x,
    seg_pair$y,
    surr,
    window_size = w,
    window_increment = w,
    ar_order = 2L
  )
  expect_equal(mv$grid$p, ref$p_value_xy)
  expect_equal(mv$grid$p_yx, ref$p_value_yx)
  expect_identical(mv$grid$window_increment, as.double(w))
})

test_that("segment rows are not crossed with increment_pct", {
  set.seed(604)
  mv <- suppressMessages(synchrony_multiverse(
    seg_pair$x,
    seg_pair$y,
    estimator = "wcc",
    sample_rate = 1,
    window_sec = 16,
    lag_sec = 2,
    increment_pct = c(0.1, 0.5),
    statistic = c("mean_abs_z", "peak"),
    surrogate_method = c("phase", "segment"),
    n_surrogates = 19
  ))
  for (stat in c("mean_abs_z", "peak")) {
    rows <- mv$grid[mv$grid$statistic == stat, ]
    expect_identical(sum(rows$surrogate_method == "phase"), 2L)
    expect_identical(sum(rows$surrogate_method == "segment"), 1L)
    seg <- rows[rows$surrogate_method == "segment", ]
    expect_true(is.na(seg$increment_pct))
    expect_identical(seg$window_increment, 20)
    phase <- rows[rows$surrogate_method == "phase", ]
    expect_setequal(phase$increment_pct, c(0.1, 0.5))
  }
})

test_that("a call with segment cells gives the aligned-design message once", {
  set.seed(605)
  n_msg <- count_conditions(
    suppressWarnings(synchrony_multiverse(
      seg_pair$x,
      seg_pair$y,
      estimator = "wcc",
      sample_rate = 1,
      window_sec = c(16, 20),
      lag_sec = 2,
      surrogate_method = "segment",
      n_surrogates = 19
    )),
    "bsync_segment_aligned"
  )
  expect_identical(n_msg, 1L)
  expect_message(
    synchrony_multiverse(
      seg_pair$x,
      seg_pair$y,
      estimator = "wcc",
      sample_rate = 1,
      window_sec = c(16, 20),
      lag_sec = 2,
      surrogate_method = "segment",
      n_surrogates = 19
    ),
    "non-overlapping windows aligned to the segments",
    class = "bsync_segment_aligned"
  )
})

test_that("segment cells with fewer than 4 segments are skipped", {
  set.seed(606)
  x <- ar1(100)
  y <- ar1(100)
  # Window 24, lag 4: segments of 32, 3 segments. Window 46, lag 4:
  # segments of 54, above floor(100 / 2) = 50.
  warnings <- list()
  mv <- withCallingHandlers(
    suppressMessages(synchrony_multiverse(
      x,
      y,
      estimator = "wcc",
      sample_rate = 1,
      window_sec = c(24, 46),
      lag_sec = 4,
      surrogate_method = "segment",
      n_surrogates = 19
    )),
    warning = function(w) {
      warnings[[length(warnings) + 1]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 2)
  for (w in warnings) {
    expect_s3_class(w, "bsync_segment_skipped")
  }
  msgs <- vapply(warnings, conditionMessage, character(1))
  expect_match(msgs[1], "`window_size` = 24 and `lag_max` = 4 samples")
  expect_match(msgs[2], "`window_size` = 46 and `lag_max` = 4 samples")
  expect_identical(nrow(mv$grid), 2L)
  expect_true(all(is.na(mv$grid$p)))
  expect_true(all(is.na(mv$grid$es)))
  expect_identical(mv$grid$window_increment, c(32, 54))
})

test_that("a segment of a single sample is skipped, not an error", {
  set.seed(609)
  # Window 1 sample: lag capped at floor(1 / 2) = 0, so the segment is 1.
  expect_warning(
    mv <- suppressMessages(synchrony_multiverse(
      seg_pair$x,
      seg_pair$y,
      estimator = "wcc",
      sample_rate = 1,
      window_sec = 1,
      lag_sec = 1,
      surrogate_method = "segment",
      n_surrogates = 19
    )),
    "segment size of 1 sample is below 2",
    class = "bsync_segment_skipped"
  )
  expect_true(is.na(mv$grid$p))
})

test_that("no aligned-design message when every segment cell is skipped", {
  set.seed(610)
  expect_no_message(
    suppressWarnings(synchrony_multiverse(
      ar1(60),
      ar1(60),
      estimator = "wcc",
      sample_rate = 1,
      window_sec = 16,
      lag_sec = 2,
      surrogate_method = "segment",
      n_surrogates = 19
    )),
    class = "bsync_segment_aligned"
  )
})

test_that("repeated user inputs keep their segment rows", {
  set.seed(611)
  mv <- suppressMessages(synchrony_multiverse(
    seg_pair$x,
    seg_pair$y,
    estimator = "wcc",
    sample_rate = 1,
    window_sec = c(16, 16),
    lag_sec = 2,
    increment_pct = c(0.1, 0.5),
    surrogate_method = c("phase", "segment"),
    n_surrogates = 19
  ))
  expect_identical(sum(mv$grid$surrogate_method == "phase"), 4L)
  expect_identical(sum(mv$grid$surrogate_method == "segment"), 2L)
})

test_that("non-segment grid columns keep their types", {
  set.seed(612)
  mv <- synchrony_multiverse(
    seg_pair$x,
    seg_pair$y,
    estimator = "wcc",
    sample_rate = 1,
    window_sec = 16,
    lag_sec = 2,
    surrogate_method = "phase",
    n_surrogates = 19
  )
  expect_type(mv$grid$window_size, "double")
  expect_type(mv$grid$window_increment, "double")
  expect_type(mv$grid$lag_max, "integer")
})

test_that("the dashboard marks no Increment tile active for segment cells", {
  skip_if_not_installed("vdiffr")
  set.seed(613)
  mv <- suppressMessages(synchrony_multiverse(
    seg_pair$x,
    seg_pair$y,
    estimator = "wcc",
    sample_rate = 1,
    window_sec = c(16, 24),
    lag_sec = 2,
    increment_pct = c(0.1, 0.5),
    surrogate_method = c("phase", "segment"),
    n_surrogates = 19
  ))
  vdiffr::expect_doppelganger("multiverse-phase-segment", function() plot(mv))
})

test_that("a skipped Granger segment cell names its window only", {
  set.seed(607)
  x <- ar1(100)
  y <- ar1(100)
  expect_warning(
    mv <- suppressMessages(synchrony_multiverse(
      x,
      y,
      estimator = "wgranger",
      sample_rate = 1,
      window_sec = 30,
      surrogate_method = "segment",
      n_surrogates = 19
    )),
    "segment cell with `window_size` = 30 samples\\.",
    class = "bsync_segment_skipped"
  )
  expect_true(is.na(mv$grid$p))
})

test_that("autotune_wcc() runs segment cells with aligned increments", {
  set.seed(608)
  # Length 100, default lags window / 2 crossed with both windows: window 10
  # gives segments of 20 (5 segments); window 24 gives segments of 34 or 48
  # (fewer than 4 segments), which are skipped in every dyad.
  dyads <- lapply(1:3, function(i) list(x = ar1(100), y = ar1(100)))
  n_msg <- 0L
  warnings <- list()
  res <- withCallingHandlers(
    autotune_wcc(
      dyads,
      sample_rate = 1,
      window_sec = c(10, 24),
      surrogate_method = "seg", # a partial name, matched as "segment"
      n_surrogates = 19,
      sig_pct = 0
    ),
    message = function(m) {
      if (inherits(m, "bsync_segment_aligned")) {
        n_msg <<- n_msg + 1L
      }
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warnings[[length(warnings) + 1]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(n_msg, 1L)
  expect_length(warnings, 1)
  expect_s3_class(warnings[[1]], "bsync_segment_skipped")
  expect_match(conditionMessage(warnings[[1]]), "dyads 1, 2, and 3")
  expect_true(any(is.na(res$dyad_multiverses[[1]]$grid$p)))
  for (mv in res$dyad_multiverses) {
    expect_true(all(mv$grid$surrogate_method == "segment"))
    expect_identical(
      mv$grid$window_increment,
      mv$grid$window_size + 2L * mv$grid$lag_max
    )
  }
})

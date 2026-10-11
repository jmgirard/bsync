# Dropping repeated cells must not change the values of the cells that stay.
# The pins below are the output of main at ab7c509 (before repeated cells were
# dropped) for the same calls and seeds, read with
# dput(as.data.frame(m$grid[, cols]), control = "digits17"). On main, the wcc
# grid had 12 cells. Cells 3, 7, and 11 (window 1 s, lag 1 s) repeated cells 1,
# 5, and 9 in samples, because the lag cap made lag 1 s into 40 samples. The
# pins keep cells 1, 2, 4, 5, 6, 8, 9, 10, and 12 of that grid.
#
# Since M019, a circular cell shifts by at least twice its own lag (D-005), so
# this two-lag grid draws two circular matrices (shifts of at least 80 and 160
# samples) where it drew one (160). That changes the random numbers of the
# circular cells and of the segment cells drawn after them. Their null_mean,
# null_sd, es, and p below (cells 4 to 9) are the output of the M019 code for
# the same call and seed. The phase cells, drawn first, and every observed and
# n_windows value still equal ab7c509.

pin_cols <- c("observed", "null_mean", "null_sd", "es", "p", "n_windows")

test_that("kept wcc cells match the pinned output (phase, circular, segment)", {
  set.seed(20261011)
  m <- suppressWarnings(suppressMessages(synchrony_multiverse(
    sim_dyad$x_A,
    sim_dyad$x_B,
    estimator = "wcc",
    sample_rate = 80,
    window_sec = c(1, 2),
    lag_sec = c(0.5, 1),
    surrogate_method = c("phase", "circular", "segment"),
    n_surrogates = 19
  )))
  pins <- data.frame(
    observed = c(
      0.091308563173046822,
      0.064729674603016554,
      0.064734583506645266,
      0.091308563173046822,
      0.064729674603016554,
      0.064734583506645266,
      0.091097658295250503,
      0.062822554004207762,
      0.065785829363197337
    ),
    null_mean = c(
      0.090381578775603266,
      0.063664778646294717,
      0.063557204178132959,
      0.090643264542892898,
      0.063357450556461539,
      0.063222060063574531,
      0.089473269149258391,
      0.063127141723722136,
      0.063932032321359938
    ),
    null_sd = c(
      0.0007423788389253847,
      0.00076342284885224038,
      0.00076011154094212578,
      0.001211120716976097,
      0.001289716248840389,
      0.00082981318727806935,
      0.0018982313185114986,
      0.0013992607515749041,
      0.0013738188358563953
    ),
    es = c(
      1.2486675923917676,
      1.3948966268468947,
      1.5489559953964065,
      0.54932478722272104,
      1.0639736048830977,
      1.822727652752872,
      0.85573824968069678,
      -0.21767759809710402,
      1.349375182122758
    ),
    p = c(0.15, 0.15, 0.10, 0.35, 0.20, 0.10, 0.25, 0.6, 0.15),
    n_windows = c(281, 136, 131, 281, 136, 131, 15, 10, 7)
  )
  expect_equal(nrow(m$grid), 9L)
  expect_equal(
    m$grid$surrogate_method,
    rep(c("phase", "circular", "segment"), each = 3)
  )
  expect_equal(m$grid$window_sec, rep(c(1, 2, 2), 3))
  expect_equal(m$grid$lag_sec, rep(c(0.5, 0.5, 1), 3))
  expect_equal(as.data.frame(m$grid[, pin_cols]), pins, tolerance = 1e-10)
})

test_that("wgranger cells without repeats match the ab7c509 output", {
  set.seed(20261012)
  g <- suppressWarnings(suppressMessages(synchrony_multiverse(
    sim_dyad$x_A,
    sim_dyad$x_B,
    estimator = "wgranger",
    sample_rate = 80,
    window_sec = c(2, 4),
    surrogate_method = "phase",
    n_surrogates = 19
  )))
  pins <- data.frame(
    observed = c(0.8841191105396351, 0.64028421616855902),
    null_mean = c(0.97480343581052709, 1.002157084236333),
    null_sd = c(0.31388660658368844, 0.44588002783766145),
    es = c(-0.28890791568933583, -0.81159245867708318),
    p = c(0.55, 0.75),
    n_windows = c(141, 66),
    observed_yx = c(2.0980982371553862, 1.5689609847998451),
    es_yx = c(3.4174710854833115, 1.0623440875692052),
    p_yx = c(0.05, 0.15)
  )
  cols <- c(pin_cols, "observed_yx", "es_yx", "p_yx")
  expect_equal(as.data.frame(g$grid[, cols]), pins, tolerance = 1e-10)
})

# A grid whose circular cells share one lag, each with at least n_surrogates
# valid shifts, draws one circular matrix with the same random numbers as
# before circular cells got their own shift range (D-005). So every cell, of
# every method, keeps its value. The pins are the output of main at c13ffca
# for the same call and seed, read with sprintf("%a", m$grid[[col]]). They are
# hexadecimal doubles because R parses a long decimal literal up to one unit
# in the last place off (0.00076342284885224038 parses one below the null_sd
# of cell 2), and this test compares with tolerance 0.
test_that("a single-lag wcc grid matches the c13ffca output exactly", {
  set.seed(20261011)
  m <- suppressWarnings(suppressMessages(synchrony_multiverse(
    sim_dyad$x_A,
    sim_dyad$x_B,
    estimator = "wcc",
    sample_rate = 80,
    window_sec = c(1, 2),
    lag_sec = 0.5,
    surrogate_method = c("phase", "circular", "segment"),
    n_surrogates = 19
  )))
  pins <- data.frame(
    observed = c(
      0x1.75fff7cac48d1p-4,
      0x1.0921fbb802c6dp-4,
      0x1.75fff7cac48d1p-4,
      0x1.0921fbb802c6dp-4,
      0x1.7522d171ecaefp-4,
      0x1.015238ee6364fp-4
    ),
    null_mean = c(
      0x1.7233f4500874p-4,
      0x1.04c55be3162ep-4,
      0x1.73465a0d07c9cp-4,
      0x1.03831a219128dp-4,
      0x1.6eb7c510fb995p-4,
      0x1.024f31c182e08p-4
    ),
    null_sd = c(
      0x1.853866ace80d5p-11,
      0x1.9040e159c8894p-11,
      0x1.3d7cef7bd28f4p-10,
      0x1.5217647086dp-10,
      0x1.24115766c62a4p-9,
      0x1.7e625c5eee1b4p-10
    ),
    es = c(
      0x1.3fa8ade8e5cfbp+0,
      0x1.6517f2019bae5p+0,
      0x1.1941193801eafp-1,
      0x1.106092fcc7b6ap+0,
      0x1.680962cc419a7p-1,
      -0x1.52b8a1bdb74f7p-3
    ),
    p = c(
      0x1.3333333333333p-3,
      0x1.3333333333333p-3,
      0x1.6666666666666p-2,
      0x1.999999999999ap-3,
      0x1p-2,
      0x1p-1
    ),
    n_windows = c(281, 136, 281, 136, 15, 10)
  )
  expect_equal(
    m$grid$surrogate_method,
    rep(c("phase", "circular", "segment"), each = 2)
  )
  expect_identical(m$grid$window_size, rep(c(80L, 160L), 3))
  expect_identical(m$grid$lag_max, rep(40L, 6))
  expect_identical(m$grid$n_surrogates, rep(19L, 6))
  expect_equal(as.data.frame(m$grid[, pin_cols]), pins, tolerance = 0)
})

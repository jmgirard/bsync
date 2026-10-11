# Dropping repeated cells must not change the values of the cells that stay.
# The pins below are the output of main at ab7c509 (before repeated cells were
# dropped) for the same calls and seeds, read with
# dput(as.data.frame(m$grid[, cols]), control = "digits17"). On main, the wcc
# grid had 12 cells. Cells 3, 7, and 11 (window 1 s, lag 1 s) repeated cells 1,
# 5, and 9 in samples, because the lag cap made lag 1 s into 40 samples. The
# pins keep cells 1, 2, 4, 5, 6, 8, 9, 10, and 12 of that grid.

pin_cols <- c("observed", "null_mean", "null_sd", "es", "p", "n_windows")

test_that("kept wcc cells match the ab7c509 output (phase, circular, segment)", {
  set.seed(20261011)
  m <- suppressWarnings(suppressMessages(synchrony_multiverse(
    sim_dyad$x_A, sim_dyad$x_B,
    estimator = "wcc", sample_rate = 80,
    window_sec = c(1, 2), lag_sec = c(0.5, 1),
    surrogate_method = c("phase", "circular", "segment"),
    n_surrogates = 19
  )))
  pins <- data.frame(
    observed = c(
      0.091308563173046822, 0.064729674603016554, 0.064734583506645266,
      0.091308563173046822, 0.064729674603016554, 0.064734583506645266,
      0.091097658295250503, 0.062822554004207762, 0.065785829363197337
    ),
    null_mean = c(
      0.090381578775603266, 0.063664778646294717, 0.063557204178132959,
      0.090164945960022377, 0.062919307059020282, 0.063251386316452857,
      0.089450197506646303, 0.063063806876179851, 0.063904795916891466
    ),
    null_sd = c(
      0.0007423788389253847, 0.00076342284885224038, 0.00076011154094212578,
      0.0014163779730158755, 0.001542345455172324, 0.00077808449811727579,
      0.0023486524150100716, 0.0014586800459551576, 0.0015712688610443693
    ),
    es = c(
      1.2486675923917676, 1.3948966268468947, 1.5489559953964065,
      0.80742374903596914, 1.1737756531295398, 1.9062160906447669,
      0.70144938351686004, -0.16539121971337753, 1.1971429542972118
    ),
    p = c(0.15, 0.15, 0.10, 0.20, 0.20, 0.05, 0.25, 0.5, 0.15),
    n_windows = c(281, 136, 131, 281, 136, 131, 15, 10, 7)
  )
  expect_equal(nrow(m$grid), 9L)
  expect_equal(m$grid$surrogate_method, rep(c("phase", "circular", "segment"), each = 3))
  expect_equal(m$grid$window_sec, rep(c(1, 2, 2), 3))
  expect_equal(m$grid$lag_sec, rep(c(0.5, 0.5, 1), 3))
  expect_equal(as.data.frame(m$grid[, pin_cols]), pins, tolerance = 1e-10)
})

test_that("wgranger cells without repeats match the ab7c509 output", {
  set.seed(20261012)
  g <- suppressWarnings(suppressMessages(synchrony_multiverse(
    sim_dyad$x_A, sim_dyad$x_B,
    estimator = "wgranger", sample_rate = 80,
    window_sec = c(2, 4), surrogate_method = "phase", n_surrogates = 19
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

# Generator for the frozen simulated-null PLV bound in test-wphase.R.
#
# Simulates the null distribution of wphase()'s mean_plv aggregate on
# independent white-noise pairs and prints its 99th percentile. The committed
# value in test-wphase.R ("white-noise PLV stays below the committed
# simulated-null q99") must equal this script's output; re-run after any
# change to the shipped PLV path.
#
# Reproducibility: master seed below; per-replicate seeds derive from it.
# Parameters match the asserting test exactly.

devtools::load_all(quiet = TRUE)

master_seed <- 20260829
n_reps <- 500
n <- 600
window_size <- 96
lag_max <- 10
window_increment <- 4
lag_increment <- 2

set.seed(master_seed)
rep_seeds <- sample.int(.Machine$integer.max, n_reps)

null_mean_plv <- vapply(rep_seeds, function(s) {
  set.seed(s)
  x <- rnorm(n)
  y <- rnorm(n)
  wphase(
    x, y,
    window_size = window_size, lag_max = lag_max,
    window_increment = window_increment, lag_increment = lag_increment
  )$aggregate[["mean_plv"]]
}, numeric(1))

q99 <- stats::quantile(null_mean_plv, 0.99, names = FALSE)
cat(sprintf(
  "n_reps = %d | mean = %.7f | q99 = %.7f\n",
  n_reps, mean(null_mean_plv), q99
))

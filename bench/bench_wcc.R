# bench/bench_wcc.R
# Before/after timing harness for the M2 prefix-sum WCC rewrite.
# Run from the package root: source("bench/bench_wcc.R")
# Requires the bench package (install.packages("bench")).

devtools::load_all(quiet = TRUE)

# ---- Configuration ----------------------------------------------------------

set.seed(42)

# Config 1: sim_dyad (n=2400, 80 Hz — the canonical reference series)
x_sim <- sim_dyad$x_A
y_sim <- sim_dyad$x_B

# Config 2: larger synthetic series (n=10000)
n_large <- 10000
x_large <- rnorm(n_large)
y_large <- rnorm(n_large)

configs <- list(
  sim_dyad_narrow  = list(x = x_sim,   y = y_sim,   window_size =  96, lag_max = 40),
  sim_dyad_wide    = list(x = x_sim,   y = y_sim,   window_size = 240, lag_max = 80),
  large_narrow     = list(x = x_large, y = y_large, window_size = 100, lag_max = 20),
  large_wide       = list(x = x_large, y = y_large, window_size = 500, lag_max = 100)
)

# ---- Benchmark --------------------------------------------------------------

results <- lapply(names(configs), function(nm) {
  cfg <- configs[[nm]]
  cat("Benchmarking:", nm, "\n")

  bm <- bench::mark(
    wcc(cfg$x, cfg$y, window_size = cfg$window_size, lag_max = cfg$lag_max),
    iterations = 5,
    check = FALSE
  )

  data.frame(
    config      = nm,
    n           = length(cfg$x),
    window_size = cfg$window_size,
    lag_max     = cfg$lag_max,
    median_ms   = as.numeric(bm$median) * 1000,
    stringsAsFactors = FALSE
  )
})

results_df <- do.call(rbind, results)
print(results_df, row.names = FALSE)

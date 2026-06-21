# data-raw/sim_dyad.R

set.seed(2026)

fs <- 80
duration <- 30
N <- fs * duration
t <- seq(0, duration, length.out = N)

freq <- 0.5

# Shifting lag: Person A leads by 0.5s, smoothly transitions to B leading by 0.5s
lag_shift_sec <- seq(0.5, -0.5, length.out = N)

# ADDED: Amplitude envelope to make the interaction non-stationary.
# Movement starts at 0, peaks in the middle of the 30s, and returns to 0.
envelope <- sin(pi * t / duration)

# Signal generation (multiplying by the envelope)
signal_A <- sin(2 * pi * freq * t) * envelope
signal_B <- sin(2 * pi * freq * (t - lag_shift_sec)) * envelope

# Simulating 3D position with highly realistic, low-level sensor jitter (sd = 0.002)
x_A <- rnorm(N, mean = 0, sd = 0.002)
y_A <- rnorm(N, mean = 0, sd = 0.002)
z_A <- signal_A + rnorm(N, mean = 0, sd = 0.002)

x_B <- rnorm(N, mean = 0, sd = 0.002)
y_B <- rnorm(N, mean = 0, sd = 0.002)
z_B <- signal_B + rnorm(N, mean = 0, sd = 0.002)

sim_dyad <- data.frame(
  time = t,
  x_A = x_A,
  y_A = y_A,
  z_A = z_A,
  x_B = x_B,
  y_B = y_B,
  z_B = z_B
)

usethis::use_data(sim_dyad, overwrite = TRUE)

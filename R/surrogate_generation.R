# =========================================================================
# === SURROGATE DATA GENERATORS ===========================================
# =========================================================================
# These functions generate null time series data to test the statistical
# significance of synchrony metrics.
# =========================================================================

# -------------------------------------------------------------------------
# --- 1. Circular Shift Method --------------------------------------------
# -------------------------------------------------------------------------

#' Generate Circular Shift Surrogates
#'
#' @param y A numeric vector containing a time series.
#' @param n_surrogates Integer specifying the number of surrogates. Default is 100.
#' @param lag_max Optional integer. If provided, ensures shifts are large enough
#'   to break local autocorrelation.
#' @return A matrix where each column is a surrogate time series.
#' @export
generate_surrogate_circular <- function(y, n_surrogates = 100, lag_max = NULL) {
  n_y <- length(y)

  if (!is.null(lag_max)) {
    min_shift <- lag_max * 2
    max_shift <- n_y - min_shift
    if (max_shift <= min_shift) {
      stop(
        "Time series is too short relative to lag_max to perform valid circular shifts."
      )
    }
    valid_shifts <- min_shift:max_shift
  } else {
    valid_shifts <- 1:(n_y - 1)
  }

  if (length(valid_shifts) < n_surrogates) {
    warning("Limited unique shifts available. Sampling with replacement.")
    shifts <- sample(valid_shifts, n_surrogates, replace = TRUE)
  } else {
    shifts <- sample(valid_shifts, n_surrogates, replace = FALSE)
  }

  # Pre-allocate and fill matrix
  surr_mat <- matrix(0, nrow = n_y, ncol = n_surrogates)
  for (j in seq_len(n_surrogates)) {
    s <- shifts[j]
    surr_mat[, j] <- c(y[(s + 1):n_y], y[1:s])
  }

  return(surr_mat)
}

# -------------------------------------------------------------------------
# --- 2. Phase Randomization Method ---------------------------------------
# -------------------------------------------------------------------------

#' Generate Phase-Randomized Surrogates (Fourier Transform)
#'
#' @param y A numeric vector containing a time series.
#' @param n_surrogates Integer specifying the number of surrogates. Default is 100.
#' @return A matrix where each column is a surrogate time series.
#' @export
generate_surrogate_phase <- function(y, n_surrogates = 100) {
  n_y <- length(y)

  # 1. Fourier transform the original signal
  y_fft <- stats::fft(y)
  amplitudes <- Mod(y_fft)

  # 2. Setup indices to ensure Hermitian symmetry (real-valued output)
  half_n <- floor(n_y / 2)

  # 3. Pre-allocate phase matrix for speed
  phase_mat <- matrix(0, nrow = n_y, ncol = n_surrogates)

  for (j in seq_len(n_surrogates)) {
    # Generate random phases uniformly from 0 to 2*pi
    random_phases <- stats::runif(half_n - 1, 0, 2 * pi)

    if (n_y %% 2 == 0) {
      # Even length
      phases <- c(0, random_phases, 0, -rev(random_phases))
    } else {
      # Odd length
      phases <- c(0, random_phases, -rev(random_phases))
    }
    phase_mat[, j] <- phases
  }

  # 4. Reconstruct complex signals with original amplitudes and random phases
  complex_mat <- amplitudes * exp(1i * phase_mat)

  # 5. Fast inverse FFT across all columns simultaneously
  surr_mat <- Re(stats::mvfft(complex_mat, inverse = TRUE)) / n_y

  return(surr_mat)
}

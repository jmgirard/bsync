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
      cli::cli_abort(
        "Time series is too short relative to {.arg lag_max} to perform valid circular shifts."
      )
    }
    valid_shifts <- min_shift:max_shift
  } else {
    valid_shifts <- seq_len(n_y - 1)
  }

  if (length(valid_shifts) < n_surrogates) {
    cli::cli_warn("Limited unique shifts available. Sampling with replacement.")
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
#' @param trim_odd Logical. If TRUE, drops the final observation if the time series length is odd.
#' @return A matrix where each column is a surrogate time series.
#' @export
generate_surrogate_phase <- function(y, n_surrogates = 100, trim_odd = FALSE) {
  n_y <- length(y)

  # 1. Handle Odd Lengths Safely
  if (n_y %% 2 != 0) {
    if (trim_odd) {
      y <- y[-n_y]
      n_y <- length(y)
      cli::cli_alert_warning(
        "Odd length detected. Trimming the final observation to {n_y}."
      )
    } else {
      cli::cli_abort(c(
        "Phase randomization requires an even number of observations.",
        "x" = "Current length is {n_y}.",
        "i" = "Set {.arg trim_odd = TRUE} to automatically drop the last observation."
      ))
    }
  }

  # 2. Fourier transform the original signal
  y_fft <- stats::fft(y)
  amplitudes <- Mod(y_fft)

  # 3. Setup indices to ensure Hermitian symmetry
  half_n <- floor(n_y / 2)

  # 4. Pre-allocate phase matrix for speed
  phase_mat <- matrix(0, nrow = n_y, ncol = n_surrogates)

  for (j in seq_len(n_surrogates)) {
    # Generate random phases uniformly from 0 to 2*pi
    random_phases <- stats::runif(half_n - 1, 0, 2 * pi)

    # Construct phases for an even-length signal
    phases <- c(0, random_phases, 0, -rev(random_phases))
    phase_mat[, j] <- phases
  }

  # 5. Reconstruct complex signals and inverse FFT
  complex_mat <- amplitudes * exp(1i * phase_mat)
  surr_mat <- Re(stats::mvfft(complex_mat, inverse = TRUE)) / n_y

  return(surr_mat)
}

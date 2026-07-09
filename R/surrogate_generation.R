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
#' @examples
#' # Build 100 circular-shift surrogates of one partner's signal
#' surr <- generate_surrogate_circular(sim_dyad$x_B, n_surrogates = 100)
#' dim(surr)
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
#' @details
#' Randomizes the phases of the Fourier transform of `y` while preserving its
#' amplitude spectrum (and therefore its autocovariance / power spectrum). The
#' DC term and -- for even-length series -- the Nyquist term are real-valued
#' components carrying a sign; they are preserved exactly (not just their
#' modulus), so the surrogate mean and variance match the original for any
#' signal, including negative-mean signals. Conjugate symmetry is enforced on
#' the remaining bins so the inverse transform is real. Both even- and
#' odd-length series are handled natively.
#'
#' @param y A numeric vector containing a time series.
#' @param n_surrogates Integer specifying the number of surrogates. Default is 100.
#' @param trim_odd Logical. Odd-length series are supported natively; when
#'   `TRUE`, the final observation is dropped instead (legacy behavior).
#'   Default is `FALSE`.
#' @return A matrix where each column is a surrogate time series.
#' @examples
#' # Build 100 phase-randomized surrogates (preserves the power spectrum)
#' surr <- generate_surrogate_phase(sim_dyad$x_B, n_surrogates = 100)
#' dim(surr)
#' @export
generate_surrogate_phase <- function(y, n_surrogates = 100, trim_odd = FALSE) {
  n_y <- length(y)

  # Legacy opt-in: drop the final observation on odd length (no longer needed
  # for correctness -- odd lengths are handled natively below -- but retained
  # for backward compatibility).
  if (n_y %% 2 != 0 && trim_odd) {
    y <- y[-n_y]
    n_y <- length(y)
    cli::cli_alert_warning(
      "Odd length detected. Trimming the final observation to {n_y}."
    )
  }

  # Fourier transform the original signal; keep amplitudes for randomized bins.
  y_fft <- stats::fft(y)
  amplitudes <- Mod(y_fft)

  n_even <- (n_y %% 2 == 0)
  half_n <- floor(n_y / 2)

  # Free positive-frequency bins whose phase we randomize (1-based indices).
  # Even n: bins 2..half_n (bin half_n + 1 is the real-valued Nyquist term).
  # Odd  n: bins 2..half_n + 1 (there is no Nyquist term).
  n_free <- if (n_even) half_n - 1L else half_n
  pos_idx <- 1L + seq_len(n_free)
  neg_idx <- n_y + 2L - pos_idx # conjugate-symmetric mirror bins

  # Build the phase-randomized complex spectrum, one column per surrogate.
  spec <- matrix(0 + 0i, nrow = n_y, ncol = n_surrogates)
  spec[1, ] <- y_fft[1] # DC term preserved exactly (keeps its sign)
  if (n_even) {
    spec[half_n + 1L, ] <- y_fft[half_n + 1L] # Nyquist term preserved exactly
  }

  amp_pos <- amplitudes[pos_idx]
  for (j in seq_len(n_surrogates)) {
    random_phases <- stats::runif(n_free, 0, 2 * pi)
    vals <- amp_pos * exp(1i * random_phases)
    spec[pos_idx, j] <- vals
    spec[neg_idx, j] <- Conj(vals)
  }

  # Inverse FFT; imaginary part is ~0 by construction.
  surr_mat <- Re(stats::mvfft(spec, inverse = TRUE)) / n_y

  return(surr_mat)
}

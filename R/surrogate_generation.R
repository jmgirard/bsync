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

# -------------------------------------------------------------------------
# --- 3. Pseudo-Dyad Method (between dyads) -------------------------------
# -------------------------------------------------------------------------
# M009. Convention: rMEA shuffle() (Kleinbub & Ramseyer, 2020) -- pairs
# series from different dyads and crops each pair to the shorter series,
# start-aligned. bsync keeps roles by default and never drops NA rows.

#' Generate Pseudo-Dyad Surrogates for One Dyad
#'
#' Builds a surrogate matrix for one target dyad whose columns are partner
#' series taken from *other* dyads in the sample. Pairing the target's `x`
#' with a partner who never interacted with that person gives a
#' pseudo-synchrony null: it keeps any co-movement that the shared task or
#' setting produces at the same moment, and removes the coupling that is
#' specific to the real interaction.
#'
#' @details
#' **Which null this tests.** Circular-shift and phase surrogates break all
#' time alignment, so a shared task structure (for example, both people move
#' when a trial starts) can look like synchrony against them. A pseudo-dyad
#' partner was recorded under the same task timeline but in a different
#' interaction, so the null keeps task-driven co-movement. A significant
#' result then means more synchrony than people in the same setting show with
#' a stranger (Kleinbub & Ramseyer, 2020).
#'
#' **Length handling.** Partners are cropped *start-aligned*: each column is
#' the first `length(y)` samples of its source, so sample `t` of the partner
#' sits at the same task time as sample `t` of the target (the rMEA
#' convention). Partners shorter than the target's `y` cannot fill the matrix
#' and are excluded with a warning. How many longer partners were cropped is
#' announced. Missing values are passed through unchanged.
#'
#' **Roles.** With `keep_roles = TRUE` (the default) every partner is the `y`
#' of another dyad. bsync's lead-lag sign depends on which series is `x`, so
#' keeping roles keeps the null comparable to the observed dyad. This differs
#' from rMEA's `shuffle()`, which mixes roles by default. Use
#' `keep_roles = FALSE` when partners are interchangeable: the `x` and the
#' `y` of every other dyad are then candidates, which doubles the pool.
#'
#' **How many surrogates.** A sample of N dyads offers at most N - 1
#' partners per target (2(N - 1) with `keep_roles = FALSE`), so the p-value
#' has a resolution of 1 / `ncol`. `n_surrogates = NULL` (the default) uses
#' every eligible partner, with no random draw. An integer draws that many
#' partners without replacement. Asking for more partners than exist is an
#' error, because repeated partners add no information to the null.
#'
#' @param dyad_list A list with one element per dyad. Each element is a data
#'   frame (first two columns are `x` and `y`) or a list (elements named `x`
#'   and `y`, else its first two elements), the form [autotune_wcc()] reads.
#' @param dyad A single integer: the index of the target dyad in `dyad_list`.
#' @param n_surrogates `NULL` (default) to use every eligible partner, or a
#'   single positive integer to draw that many without replacement.
#' @param keep_roles Logical. `TRUE` (default) draws partners only from the
#'   `y` series of other dyads. `FALSE` also allows their `x` series.
#' @return A numeric matrix with `length(y)` rows and one column per partner,
#'   ready to pass as `y_surrogates` to [wcc_surrogate()] and the other
#'   surrogate wrappers. The attribute `"sources"` is a data frame with
#'   columns `dyad` and `role` giving each column's source.
#' @references Kleinbub, J. R., & Ramseyer, F. T. (2020). rMEA: An R package
#'   to assess nonverbal synchronization in motion energy analysis
#'   time-series. *Psychotherapy Research*. \doi{10.1080/10503307.2020.1844334}
#' @seealso [generate_pseudo_dyads()] for the sample-wide set of
#'   pseudo-dyads; [generate_surrogate_circular()] and
#'   [generate_surrogate_phase()] for within-dyad nulls; [wcc_surrogate()],
#'   [wdtw_surrogate()], [wgranger_surrogate()], [wphase_surrogate()].
#' @examples
#' # Three "dyads" built from sim_dyad's axes (a stand-in for a real sample)
#' dyads <- list(
#'   list(x = sim_dyad$x_A, y = sim_dyad$x_B),
#'   list(x = sim_dyad$y_A, y = sim_dyad$y_B),
#'   list(x = sim_dyad$z_A, y = sim_dyad$z_B)
#' )
#' y_pseudo <- generate_surrogate_pseudo(dyads, dyad = 3)
#' dim(y_pseudo)
#' attr(y_pseudo, "sources")
#'
#' \donttest{
#' # Test dyad 3 against partners from the other dyads
#' wcc_surrogate(
#'   x = sim_dyad$z_A, y = sim_dyad$z_B, y_surrogates = y_pseudo,
#'   window_size = 96, lag_max = 10
#' )
#' }
#' @export
generate_surrogate_pseudo <- function(dyad_list,
                                      dyad,
                                      n_surrogates = NULL,
                                      keep_roles = TRUE) {
  dyads <- .pseudo_extract_dyads(dyad_list)
  if (!.is_count(dyad) || dyad > length(dyads)) {
    cli::cli_abort(
      "{.arg dyad} must be a single index into {.arg dyad_list} \\
      (1 to {length(dyads)})."
    )
  }
  if (!is.null(n_surrogates) && !.is_count(n_surrogates)) {
    cli::cli_abort(
      "{.arg n_surrogates} must be NULL or a single positive integer."
    )
  }
  .check_keep_roles(keep_roles)

  y <- dyads[[dyad]]$y
  n_y <- length(y)

  # Candidate partner series from every other dyad.
  roles <- if (keep_roles) "y" else c("x", "y")
  pool <- expand.grid(
    role = roles, dyad = setdiff(seq_along(dyads), dyad),
    stringsAsFactors = FALSE
  )[, c("dyad", "role")]
  pool_len <- mapply(
    function(d, r) length(dyads[[d]][[r]]), pool$dyad, pool$role
  )

  eligible <- pool_len >= n_y
  n_excluded <- sum(!eligible)
  if (n_excluded > 0) {
    cli::cli_warn(
      "Excluded {n_excluded} partner series shorter than the target's \\
      {.arg y} ({n_y} samples)."
    )
  }
  if (!any(eligible)) {
    cli::cli_abort(
      "No partner series is at least as long as the target's {.arg y} \\
      ({n_y} samples)."
    )
  }
  pool <- pool[eligible, , drop = FALSE]
  n_cropped <- sum(pool_len[eligible] > n_y)

  n_eligible <- nrow(pool)
  if (is.null(n_surrogates)) {
    pick <- seq_len(n_eligible)
  } else {
    if (n_surrogates > n_eligible) {
      cli::cli_abort(
        "{.arg n_surrogates} ({n_surrogates}) exceeds the {n_eligible} \\
        eligible partner series."
      )
    }
    pick <- sample.int(n_eligible, n_surrogates)
  }
  pool <- pool[pick, , drop = FALSE]

  if (n_cropped > 0) {
    cli::cli_inform(
      "Cropped {n_cropped} partner series to their first {n_y} samples \\
      (start-aligned)."
    )
  }

  surr_mat <- matrix(0, nrow = n_y, ncol = nrow(pool))
  for (j in seq_len(nrow(pool))) {
    surr_mat[, j] <- dyads[[pool$dyad[j]]][[pool$role[j]]][seq_len(n_y)]
  }
  rownames(pool) <- NULL
  attr(surr_mat, "sources") <- pool

  surr_mat
}

# --- Pseudo-dyad helpers (internal) --------------------------------------

# Extract and type-check every dyad in a dyad_list (>= 2 dyads).
.pseudo_extract_dyads <- function(dyad_list) {
  if (!is.list(dyad_list) || is.data.frame(dyad_list) ||
    length(dyad_list) < 2) {
    cli::cli_abort(
      "{.arg dyad_list} must contain at least two dyads (a list with one \\
      element per dyad)."
    )
  }
  lapply(seq_along(dyad_list), function(i) {
    xy <- .extract_xy(dyad_list[[i]])
    if (!is.numeric(xy$x) || !is.numeric(xy$y)) {
      cli::cli_abort(
        "Both series of dyad {i} in {.arg dyad_list} must be numeric."
      )
    }
    list(x = as.double(xy$x), y = as.double(xy$y))
  })
}

# TRUE for a single, non-missing, positive whole number.
.is_count <- function(v) {
  rlang::is_integerish(v, n = 1, finite = TRUE) && v >= 1
}

.check_keep_roles <- function(keep_roles) {
  if (!rlang::is_bool(keep_roles)) {
    cli::cli_abort("{.arg keep_roles} must be TRUE or FALSE.")
  }
}

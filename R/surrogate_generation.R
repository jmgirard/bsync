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
#' @seealso \code{\link{generate_surrogate_phase}},
#'   \code{\link{generate_surrogate_iaaft}}, and
#'   \code{\link{generate_surrogate_segment}} for the other within-dyad nulls;
#'   \code{\link{generate_surrogate_pseudo}} for the between-dyad null.
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
#' @seealso \code{\link{generate_surrogate_iaaft}}, which also keeps the
#'   value distribution, \code{\link{generate_surrogate_circular}}, and
#'   \code{\link{generate_surrogate_segment}} for the other within-dyad
#'   nulls; \code{\link{generate_surrogate_pseudo}}
#'   for the between-dyad null.
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
# --- 3. IAAFT Method -----------------------------------------------------
# -------------------------------------------------------------------------
# M014. Source: Schreiber & Schmitz (1996), p. 2 (cairn/references/
# schreiber1996.md).

#' Generate IAAFT Surrogates
#'
#' Builds surrogates with the iterative amplitude-adjusted Fourier transform
#' (IAAFT) of Schreiber and Schmitz (1996). Each surrogate is a new series
#' that matches both the power spectrum and the value distribution of `y`.
#' One of the two matches exactly and the other closely, as `match` selects.
#'
#' @details
#' **Which null this tests.** Phase randomization ([generate_surrogate_phase()])
#' keeps the power spectrum but makes the values Gaussian, so a skewed or
#' bounded signal gets surrogates with a different value distribution. IAAFT
#' keeps the spectrum and the value distribution. The null is a Gaussian
#' linear process seen through a fixed, monotone transform, for example a
#' skewed movement-energy signal. Against this null, a significant result
#' cannot come from the shape of the value distribution or from the
#' autocorrelation of `y` alone.
#'
#' **Algorithm.** Each surrogate starts from a random shuffle of `y`. Each
#' iteration takes the Fourier transform, replaces its amplitudes with those
#' of `y` while it keeps the phases, and transforms back. This is the
#' spectrum-adjusted series. The iteration then rank-orders it, so that the
#' series takes exactly the values of `y`. A surrogate stops when the
#' rank-ordering no longer changes it, which is where the iteration ends
#' (Schreiber & Schmitz, 1996, p. 2). Apart from a cyclic shift of `y`, a
#' finite series in general cannot match the spectrum and the values both
#' exactly (p. 2). A very short series has few orderings, and its surrogates
#' can be cyclic shifts or reversals of `y`, which break little of its
#' structure.
#'
#' **`match`.** `"spectrum"` (the default) returns the last spectrum-adjusted
#' series. Its Fourier amplitudes equal those of `y` exactly, and its values
#' are close to the values of `y`: closer than the values of a phase
#' surrogate, in the package tests on a skewed autoregressive series.
#' `"values"` returns the rank-ordered series, as in the paper. It holds
#' exactly the values of `y`, and its spectrum is close to that of `y`. The
#' rank step lowers the autocorrelation slightly, and a synchrony statistic
#' that grows with autocorrelation, such as the mean absolute Fisher z of
#' [wcc()], then reads high against these surrogates. So a test against
#' `"values"` surrogates can reject a true null more often than its nominal
#' level, and the spectrum form is the default. In a size check run on
#' 2026-10-06 (independent pairs of AR(1) series x_n = 0.7 x_(n-1) + e_n
#' observed as x_n^3, length 512, 19 surrogates per pair, [wcc_surrogate()]
#' with `window_size = 32`, `lag_max = 4`, `window_increment = 16`), the
#' spectrum form rejected 17 of 400 pairs (4.25%) at p <= .05. On the first
#' 200 of those pairs, the values form rejected 21 (10.5%) and the spectrum
#' form 9 (4.5%). The package's own tests repeat this check.
#'
#' **`max_iter`.** The spectral error falls about as 1 / i over the first
#' iterations (p. 3), and the paper's Fig. 2 runs to 1000 iterations (p. 2).
#' The default `1000` is a cap on run time, not a target. A surrogate that
#' reaches `max_iter` without a fixed point is still returned, and the call
#' warns with the number of such surrogates. Long series need more
#' iterations: on one AR(1)-cube series of 10000 samples (run 2026-10-06),
#' 4 of 5 surrogates had not reached the fixed point at 1000 iterations.
#' With the default `match`, such surrogates still match the spectrum
#' exactly.
#'
#' **Missing values.** The Fourier transform cannot take `NA`, so `y` must
#' be complete. Fill gaps first, for example with [impute_ts_gaps()].
#'
#' @param y A numeric vector with at least 3 finite values and no `NA`.
#' @param n_surrogates A single positive integer: the number of surrogates.
#'   Default is `100`.
#' @param max_iter A single positive integer: the maximum number of
#'   iterations per surrogate. Default is `1000`.
#' @param match Which property matches exactly: `"spectrum"` (default), the
#'   Fourier amplitudes of `y`, or `"values"`, the values of `y`. The other
#'   property matches closely. See Details.
#' @return A numeric matrix with `length(y)` rows and `n_surrogates`
#'   columns, ready to pass as `y_surrogates` to [wcc_surrogate()] and the
#'   other surrogate wrappers.
#' @references Schreiber, T., & Schmitz, A. (1996). Improved surrogate data
#'   for nonlinearity tests. *Physical Review Letters*, 77(4), 635-638.
#'   \doi{10.1103/PhysRevLett.77.635}
#' @seealso [generate_surrogate_phase()], [generate_surrogate_circular()],
#'   and [generate_surrogate_segment()] for the other within-dyad nulls;
#'   [generate_surrogate_pseudo()] for the
#'   between-dyad null; [wcc_surrogate()], [wdtw_surrogate()],
#'   [wgranger_surrogate()], [wphase_surrogate()].
#' @examples
#' # Build 100 IAAFT surrogates: exact spectrum, close values
#' surr <- generate_surrogate_iaaft(sim_dyad$x_B, n_surrogates = 100)
#' dim(surr)
#' all.equal(Mod(fft(surr[, 1])), Mod(fft(sim_dyad$x_B)))
#'
#' # The "values" form: exact values, close spectrum
#' surr_v <- generate_surrogate_iaaft(sim_dyad$x_B, 10, match = "values")
#' all.equal(sort(surr_v[, 1]), sort(sim_dyad$x_B))
#' @md
#' @export
generate_surrogate_iaaft <- function(
  y,
  n_surrogates = 100,
  max_iter = 1000,
  match = c("spectrum", "values")
) {
  # arg_match() takes the first element of any reordering of the choices, so
  # require one value unless the default vector is passed unchanged.
  if (length(match) > 1 && !identical(match, c("spectrum", "values"))) {
    cli::cli_abort(
      "{.arg match} must be a single value: {.val spectrum} or {.val values}."
    )
  }
  match <- rlang::arg_match(match)
  if (!is.numeric(y) || !is.null(dim(y))) {
    cli::cli_abort("{.arg y} must be a numeric vector.")
  }
  if (anyNA(y)) {
    cli::cli_abort(c(
      "{.arg y} must not contain missing values.",
      "i" = "Fill gaps first, for example with {.fn impute_ts_gaps}."
    ))
  }
  if (!all(is.finite(y))) {
    cli::cli_abort("{.arg y} must contain only finite values.")
  }
  if (length(y) < 3) {
    cli::cli_abort("{.arg y} must have at least 3 values.")
  }
  if (!.is_count(n_surrogates)) {
    cli::cli_abort("{.arg n_surrogates} must be a single positive integer.")
  }
  if (!.is_count(max_iter)) {
    cli::cli_abort("{.arg max_iter} must be a single positive integer.")
  }

  y <- as.double(y)
  n_y <- length(y)
  y_sorted <- sort(y)
  target_amp <- Mod(stats::fft(y))
  if (!all(is.finite(target_amp))) {
    cli::cli_abort(c(
      "The Fourier transform of {.arg y} overflows.",
      "i" = "Rescale {.arg y}, for example to standard units, first."
    ))
  }

  # Start: one random shuffle per surrogate, column 1 first.
  surr_mat <- matrix(0, nrow = n_y, ncol = n_surrogates)
  for (j in seq_len(n_surrogates)) {
    surr_mat[, j] <- y[sample.int(n_y)]
  }

  # Iterate the columns that have not yet reached a fixed point. `adj_mat`
  # keeps each column's last spectrum-adjusted series (before the rank step).
  adj_mat <- matrix(0, nrow = n_y, ncol = n_surrogates)
  active <- seq_len(n_surrogates)
  # A double counter, so a max_iter above the integer range cannot overflow.
  iter <- 0
  while (length(active) > 0 && iter < max_iter) {
    iter <- iter + 1
    current <- surr_mat[, active, drop = FALSE]
    spec <- stats::mvfft(current)
    spec <- target_amp * exp(1i * Arg(spec))
    adjusted <- Re(stats::mvfft(spec, inverse = TRUE)) / n_y
    # Rank within each column in one stable order() call; stability breaks
    # ties by position, as rank(ties.method = "first") does.
    n_active <- length(active)
    ord <- order(rep(seq_len(n_active), each = n_y), adjusted)
    ranked <- integer(n_y * n_active)
    ranked[ord] <- rep.int(seq_len(n_y), n_active)
    updated <- matrix(y_sorted[ranked], nrow = n_y)
    changed <- colSums(updated != current) > 0
    adj_mat[, active] <- adjusted
    surr_mat[, active] <- updated
    active <- active[changed]
  }

  if (length(active) > 0) {
    hint <- if (match == "spectrum") {
      "values closer to the values of {.arg y}"
    } else {
      "a spectrum closer to the spectrum of {.arg y}"
    }
    n_unconverged <- length(active)
    cli::cli_warn(
      c(
        "{n_unconverged} of {n_surrogates} surrogate{?s} did not converge \\
        within {.arg max_iter} = {max_iter} iterations.",
        "i" = paste0(
          "They are returned as they are. Raise {.arg max_iter} for ",
          hint,
          "."
        )
      ),
      class = "bsync_iaaft_unconverged",
      n_unconverged = n_unconverged,
      n_surrogates = n_surrogates,
      max_iter = max_iter
    )
  }

  if (match == "spectrum") adj_mat else surr_mat
}

# -------------------------------------------------------------------------
# --- 4. Segment-Shuffling Method -----------------------------------------
# -------------------------------------------------------------------------
# M015. Convention: SUSY 0.1.0 susy() (Tschacher & Meier, 2020) -- k =
# floor(n / segment_size) segments cut from sample 1 (cairn/references/
# tschacher2020.md). bsync reorders the segments into a full-length series
# and keeps the tail in place, where SUSY pairs segments and drops the tail.

#' Generate Segment-Shuffling Surrogates
#'
#' Cuts `y` into segments of `segment_size` samples and builds each
#' surrogate by putting the segments in a new order. Each segment keeps its
#' own values, autocorrelation, and any trend inside it. Only the time at
#' which each segment occurred changes.
#'
#' @details
#' **Which null this tests.** Each window of `x` has one segment of `y` that
#' occurred alongside it. The test asks whether that segment resembles the
#' window more than the other segments of `y` do. Under the null, the time
#' at which each segment of `y` occurred carries no information about `x`,
#' so any order of the segments is as good as the real one. The test is
#' exact if the segments of `y` are exchangeable, which is close to true
#' when a segment is long relative to the memory of `y`, and if every window
#' lies inside one segment (see "Windows and segments").
#'
#' **Segments and orders.** The series has k = `floor(length(y) /
#' segment_size)` segments, cut from sample 1. Segment i covers samples
#' (i - 1) * `segment_size` + 1 to i * `segment_size`. Each surrogate takes
#' an order drawn uniformly from all k! - 1 orders other than the original,
#' and no order is drawn twice. So with 3 or more segments, a surrogate can
#' leave a segment in place: on average (k! - k) / (k! - 1) of the k
#' segments, which is close to one from 4 segments on. This rule makes the
#' add-one p-value of the surrogate wrappers exact. From 4 segments on,
#' orders that move every segment do not form a group with the original
#' order, and a test against them rejects a true null too often. With up
#' to 8 segments and an `n_surrogates` of at least k! - 1, the call returns
#' all k! - 1 orders and says so. With 3 segments or fewer, k! - 1 is
#' below 19, so no test can reach p <= .05, and the call warns.
#'
#' **The tail.** The last `length(y) - k * segment_size` samples do not fill
#' a segment. They stay in place at the end of every surrogate, and the call
#' says how many there are. Missing values move with their segments.
#'
#' **Windows and segments.** A window that crosses a segment boundary is
#' continuous in `y` but not in a surrogate, and the test then rejects a true
#' null too often. Every window and every lagged window lies inside one
#' segment under these settings:
#'
#' * For [wcc()], [wdtw()], and [wphase()]: `segment_size >= window_size +
#'   2 * lag_max` and `window_increment = segment_size`.
#' * For [wgranger()]: `segment_size = window_size = window_increment`.
#'
#' The window grid then has k - 1 windows, or k windows when the tail is at
#' least `window_size + 2 * lag_max - 1` samples (`window_size - 1` for
#' [wgranger()]). With `segment_size = window_size + 2 * lag_max`, that is a
#' tail of `segment_size - 1` samples. With k - 1 windows, the last segment
#' serves only as a source for the surrogates. [wphase()] takes the Hilbert
#' transform of the whole series, so a segment test of phase synchrony is
#' approximate even with aligned windows. [synchrony_multiverse()] and
#' [autotune_wcc()] use these settings for `surrogate_method = "segment"`.
#'
#' **SUSY.** The segment cut follows the `susy()` function of the SUSY
#' package (Tschacher & Meier, 2020). SUSY correlates each segment of one
#' series with each segment of the other and compares the real pairs with
#' the pairs of different segments. bsync instead builds full-length
#' reordered series, so that every surrogate wrapper can use them. SUSY's
#' pseudo pairs never match a segment with itself, but bsync's orders can
#' leave a segment in place, which keeps the p-value exact. SUSY drops the
#' tail, and bsync keeps it. bsync does not reproduce SUSY's effect size or
#' its numbers. SUSY's rule `segment >= 2 * maxlag` has the analogue
#' `segment_size >= window_size + 2 * lag_max`.
#'
#' @param y A numeric vector containing a time series. `NA` values are
#'   allowed.
#' @param segment_size A single whole number from 2 to
#'   `floor(length(y) / 2)`: the segment length in samples.
#' @param n_surrogates A single positive integer: the number of surrogates.
#'   Default is `100`. Above 8 segments, it must be below k! - 1.
#' @return A numeric matrix with `length(y)` rows and one column per
#'   surrogate, ready to pass as `y_surrogates` to [wcc_surrogate()] and the
#'   other surrogate wrappers. It has `n_surrogates` columns, or k! - 1
#'   columns if that is fewer.
#' @references Tschacher, W., & Meier, D. (2020). Physiological synchrony in
#'   psychotherapy sessions. *Psychotherapy Research*, 30(5), 558-573.
#'   \doi{10.1080/10503307.2019.1612114}
#' @seealso [generate_surrogate_circular()], [generate_surrogate_phase()],
#'   and [generate_surrogate_iaaft()] for the other within-dyad nulls;
#'   [generate_surrogate_pseudo()] for the between-dyad null;
#'   [wcc_surrogate()], [wdtw_surrogate()], [wgranger_surrogate()],
#'   [wphase_surrogate()].
#' @examples
#' # Windows of 96 samples with lags up to 10 need segments of at least
#' # 96 + 2 * 10 = 116 samples: 20 segments of sim_dyad's 2400 samples.
#' surr <- generate_surrogate_segment(
#'   sim_dyad$x_B,
#'   segment_size = 116,
#'   n_surrogates = 19
#' )
#' dim(surr)
#'
#' \donttest{
#' # Step the windows one segment at a time, so each lies inside a segment
#' wcc_surrogate(
#'   x = sim_dyad$x_A, y = sim_dyad$x_B, y_surrogates = surr,
#'   window_size = 96, lag_max = 10, window_increment = 116
#' )
#' }
#' @md
#' @export
generate_surrogate_segment <- function(y, segment_size, n_surrogates = 100) {
  if (!is.numeric(y) || !is.null(dim(y))) {
    cli::cli_abort("{.arg y} must be a numeric vector.")
  }
  if (!.is_count(n_surrogates)) {
    cli::cli_abort("{.arg n_surrogates} must be a single positive integer.")
  }
  if (!.is_count(segment_size) || segment_size < 2) {
    cli::cli_abort(
      "{.arg segment_size} must be a single whole number of at least 2."
    )
  }
  n_y <- length(y)
  half <- floor(n_y / 2)
  if (segment_size > half) {
    cli::cli_abort(c(
      "{.arg segment_size} ({segment_size}) must be at most \\
      {.code floor(length(y) / 2)} ({half}).",
      "i" = "The series needs at least two segments."
    ))
  }

  y <- as.double(y)
  segment_size <- as.integer(segment_size)
  k <- n_y %/% segment_size
  # Every order except the original (D-004): with the original order they
  # form a group, which makes the add-one p-value exact.
  n_orders <- factorial(k) - 1

  if (n_orders < 19) {
    cli::cli_warn(
      c(
        "{k} segments have only {n_orders} order{?s} other than the \\
        original, so no test on these surrogates can reach p <= .05.",
        "i" = "Use a smaller {.arg segment_size}, so that the series has at \\
        least 4 segments."
      ),
      class = "bsync_segment_few_orders"
    )
  }

  if (k <= 8) {
    # Few segments: list every order except the original, then draw distinct
    # rows from the list.
    orders <- .segment_orders(k)
    if (n_surrogates >= n_orders) {
      cli::cli_inform(
        "Returning all {n_orders} possible segment order{?s}, because \\
        {.arg n_surrogates} ({n_surrogates}) is not below that number.",
        class = "bsync_segment_all_orders"
      )
    } else {
      orders <- orders[sample.int(n_orders, n_surrogates), , drop = FALSE]
    }
  } else {
    # Many segments: draw orders and reject the original order or a repeat.
    if (n_surrogates >= n_orders) {
      cli::cli_abort(
        "{.arg n_surrogates} ({n_surrogates}) must be below the \\
        {n_orders} possible orders of {k} segments."
      )
    }
    orders <- matrix(0L, nrow = n_surrogates, ncol = k)
    seen <- new.env(hash = TRUE, size = n_surrogates)
    n_found <- 0L
    while (n_found < n_surrogates) {
      p <- sample.int(k)
      if (all(p == seq_len(k))) {
        next
      }
      key <- paste(p, collapse = ",")
      if (!is.null(seen[[key]])) {
        next
      }
      seen[[key]] <- TRUE
      n_found <- n_found + 1L
      orders[n_found, ] <- p
    }
  }

  n_tail <- n_y - k * segment_size
  if (n_tail > 0) {
    cli::cli_inform(
      "The last {n_tail} sample{?s} of {.arg y} do{?es/} not fill a whole \\
      segment and stay in place.",
      class = "bsync_segment_tail"
    )
  }

  # Row index of each surrogate sample in y: segment orders[j, i] of y goes
  # to position i, and the tail rows map to themselves.
  offsets <- seq_len(segment_size)
  tail_rows <- seq_len(n_tail) + k * segment_size
  surr_mat <- matrix(0, nrow = n_y, ncol = nrow(orders))
  for (j in seq_len(nrow(orders))) {
    rows <- c(
      rep((orders[j, ] - 1L) * segment_size, each = segment_size) + offsets,
      tail_rows
    )
    surr_mat[, j] <- y[rows]
  }
  surr_mat
}

# Every order of 1..k except the original, one per row, in lexicographic
# order. Built one position at a time: each row extends by every unused
# value. The original order is the first row, so it is dropped at the end.
.segment_orders <- function(k) {
  orders <- matrix(integer(0), nrow = 1, ncol = 0)
  for (pos in seq_len(k)) {
    n_r <- nrow(orders)
    ok <- matrix(TRUE, nrow = n_r, ncol = k)
    for (c in seq_len(ncol(orders))) {
      ok[cbind(seq_len(n_r), orders[, c])] <- FALSE
    }
    idx <- which(ok, arr.ind = TRUE)
    idx <- idx[order(idx[, "row"], idx[, "col"]), , drop = FALSE]
    orders <- cbind(orders[idx[, "row"], , drop = FALSE], idx[, "col"])
  }
  unname(orders[-1, , drop = FALSE])
}

# -------------------------------------------------------------------------
# --- 5. Pseudo-Dyad Method (between dyads) -------------------------------
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
#' announced. Missing values are passed through unchanged. Start alignment
#' assumes that every dyad was sampled at the same rate and that sample 1 of
#' every series is the same task onset. Resample and trim first if not; this
#' function cannot check it.
#'
#' **Roles.** With `keep_roles = TRUE` (the default) every partner is the `y`
#' of another dyad. bsync's lead-lag sign depends on which series is `x`, so
#' keeping roles keeps the null comparable to the observed dyad. This differs
#' from rMEA's `shuffle()`, which mixes roles by default. Use
#' `keep_roles = FALSE` when partners are interchangeable: the `x` and the
#' `y` of every other dyad are then candidates, which doubles the pool.
#'
#' **How many surrogates.** A sample of N dyads offers at most N - 1
#' partners per target (2(N - 1) with `keep_roles = FALSE`). The wrappers
#' report the add-one p-value (b + 1) / (n + 1), where n is the number of
#' partner columns and b counts the partners that score at least as high as
#' the real partner. The smallest possible p-value is 1 / (n + 1). With N
#' dyads and `keep_roles = TRUE`, n is at most N - 1, so the single-dyad
#' p-value is never below 1 / N, and `p <= .05` needs at least 20 dyads.
#' With `keep_roles = FALSE`, n is at most 2(N - 1), so the floor is
#' 1 / (2N - 1) and `p <= .05` needs at least 11 dyads. When every partner is used with
#' `keep_roles = TRUE`, the p-value is the real partner's rank among the N
#' series divided by N. If the real partner is no different from a stranger,
#' each rank has chance 1 / N, so the test keeps its nominal size. With too
#' few dyads to reach `p <= .05`, prefer the sample-level comparison of
#' [generate_pseudo_dyads()]. `n_surrogates = NULL` (the default) uses
#' every eligible partner, with no random draw. An integer draws that many
#' partners without replacement. Asking for more partners than exist is an
#' error, because repeated partners add no information to the null.
#'
#' @param dyad_list A list with one element per dyad, the form
#'   [autotune_wcc()] also reads. Each element is a data frame or a list. If
#'   its names contain `x` and `y` exactly once each, those two elements are
#'   read by name, in any position. If its names contain neither `x` nor `y`,
#'   or it has no names, its first two elements are read as `x` and `y`. Any
#'   other use of the names `x` and `y` is an error. Names are never matched
#'   partially. If `dyad_list` itself is named, the names identify the dyads
#'   in the output, so every name must be non-empty, not `NA`, and unique.
#' @param dyad A single integer: the index of the target dyad in `dyad_list`.
#' @param n_surrogates `NULL` (default) to use every eligible partner, or a
#'   single positive integer to draw that many without replacement.
#' @param keep_roles Logical. `TRUE` (default) draws partners only from the
#'   `y` series of other dyads. `FALSE` also allows their `x` series.
#' @return A numeric matrix with `length(y)` rows and one column per partner,
#'   ready to pass as `y_surrogates` to [wcc_surrogate()] and the other
#'   surrogate wrappers. The attribute `"sources"` is a data frame with one
#'   row per column and the columns `dyad` (the source dyad's index in
#'   `dyad_list`), `role` (`"x"` or `"y"`), and `name` (the source dyad's name
#'   in `dyad_list`, or `NA` if `dyad_list` has no names). If `dyad_list` is
#'   named, the matrix column names are `"<name>:<role>"`. Subsetting the
#'   matrix drops the attribute, so read it before you subset.
#' @references Kleinbub, J. R., & Ramseyer, F. T. (2020). rMEA: An R package
#'   to assess nonverbal synchronization in motion energy analysis
#'   time-series. *Psychotherapy Research*. \doi{10.1080/10503307.2020.1844334}
#' @seealso [generate_pseudo_dyads()] for the sample-wide set of
#'   pseudo-dyads; [generate_surrogate_circular()],
#'   [generate_surrogate_phase()], [generate_surrogate_iaaft()], and
#'   [generate_surrogate_segment()] for within-dyad nulls; [wcc_surrogate()],
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
#' @md
#' @export
generate_surrogate_pseudo <- function(
  dyad_list,
  dyad,
  n_surrogates = NULL,
  keep_roles = TRUE
) {
  dyads <- .pseudo_extract_dyads(dyad_list)
  dyad_names <- .pseudo_dyad_names(dyad_list)
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
    role = roles,
    dyad = setdiff(seq_along(dyads), dyad),
    stringsAsFactors = FALSE
  )[, c("dyad", "role")]
  pool_len <- mapply(
    function(d, r) length(dyads[[d]][[r]]),
    pool$dyad,
    pool$role
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
  pool_len <- pool_len[eligible]

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

  # Count crops among the partners actually returned.
  n_cropped <- sum(pool_len[pick] > n_y)
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
  pool$name <- if (is.null(dyad_names)) {
    rep(NA_character_, nrow(pool))
  } else {
    dyad_names[pool$dyad]
  }
  if (!is.null(dyad_names)) {
    colnames(surr_mat) <- paste0(pool$name, ":", pool$role)
  }
  attr(surr_mat, "sources") <- pool

  surr_mat
}

#' Generate the Sample-Wide Set of Pseudo-Dyads
#'
#' Pairs series from *different* dyads to build pseudo-dyads: pairs of people
#' who never interacted. Computing the same synchrony statistic on the real
#' dyads and on these pseudo-dyads gives the classic pseudo-synchrony
#' comparison: real dyads should be more synchronous than pseudo-dyads
#' recorded under the same task (Kleinbub & Ramseyer, 2020).
#'
#' @details
#' **Pairings.** With `keep_roles = TRUE` (the default), each pseudo-dyad is
#' the `x` of dyad i with the `y` of dyad j, for every i != j: N(N - 1)
#' ordered pairs for N dyads. bsync's lead-lag sign depends on which series is
#' `x`, so keeping roles keeps pseudo-dyads comparable to real ones. With
#' `keep_roles = FALSE`, any two series from different dyads form a pair, in
#' either role: 2N(N - 1) unordered pairs. This is the pair set of rMEA's
#' `shuffle()`, which mixes roles by default, with the same assignment of
#' series to `x` and `y`: a pair of an `x` series and a `y` series keeps the
#' `x` series as `x`.
#'
#' **Length handling.** Both series of a pseudo-dyad are cropped to the
#' shorter one, start-aligned (the rMEA convention), so sample `t` of each
#' series sits at the same task time. The number of cropped pseudo-dyads is
#' announced. Missing values are passed through unchanged. Start alignment
#' assumes that every dyad was sampled at the same rate and that sample 1 of
#' every series is the same task onset. Resample and trim first if not; this
#' function cannot check it.
#'
#' **How many.** `n_pairs = NULL` (the default) returns every pairing in a
#' fixed order, with no random draw. An integer draws that many pairings
#' without replacement. Asking for more than exist is an error.
#'
#' @inheritParams generate_surrogate_pseudo
#' @param n_pairs `NULL` (default) to return every pairing, or a single
#'   positive integer to draw that many without replacement.
#' @param keep_roles Logical. `TRUE` (default) pairs the `x` of one dyad with
#'   the `y` of another. `FALSE` pairs any two series from different dyads.
#' @return A list of pseudo-dyads. Each element is `list(x = , y = )`, the
#'   dyad form that [autotune_wcc()] and [generate_surrogate_pseudo()] read.
#'   The attribute `"sources"` is a data frame with one row per pseudo-dyad
#'   and the columns `x_dyad`, `x_role`, `x_name`, `y_dyad`, `y_role`, and
#'   `y_name`: the index, role, and name in `dyad_list` of each series. The
#'   name columns are `NA` if `dyad_list` has no names. If `dyad_list` is
#'   named, the list elements are named
#'   `"<x_name>:<x_role>|<y_name>:<y_role>"`. These element names are labels:
#'   if dyad names contain `:` or `|`, two labels can be equal, and the list
#'   cannot then be passed back as a named `dyad_list`. Subsetting the list
#'   drops the attribute, so read it before you subset.
#' @references Kleinbub, J. R., & Ramseyer, F. T. (2020). rMEA: An R package
#'   to assess nonverbal synchronization in motion energy analysis
#'   time-series. *Psychotherapy Research*. \doi{10.1080/10503307.2020.1844334}
#' @seealso [generate_surrogate_pseudo()] for a per-dyad surrogate matrix,
#'   which [wcc_surrogate()], [wdtw_surrogate()], [wgranger_surrogate()], and
#'   [wphase_surrogate()] accept; [generate_surrogate_circular()],
#'   [generate_surrogate_phase()], [generate_surrogate_iaaft()], and
#'   [generate_surrogate_segment()] for within-dyad nulls; [wcc()], [wdtw()],
#'   [wgranger()], and [wphase()] to compute the statistic on each
#'   pseudo-dyad.
#' @examples
#' # Three "dyads" built from sim_dyad's axes (a stand-in for a real sample)
#' dyads <- list(
#'   list(x = sim_dyad$x_A, y = sim_dyad$x_B),
#'   list(x = sim_dyad$y_A, y = sim_dyad$y_B),
#'   list(x = sim_dyad$z_A, y = sim_dyad$z_B)
#' )
#' pseudo <- generate_pseudo_dyads(dyads)
#' length(pseudo) # 3 * 2 = 6 role-keeping pairings
#' attr(pseudo, "sources")
#'
#' # Mean |Fisher z| on each pseudo-dyad: the pseudo-synchrony baseline
#' pseudo_z <- vapply(pseudo, function(d) {
#'   wcc(d$x, d$y, window_size = 96, lag_max = 10, window_increment = 48)$
#'     aggregate[[1]]
#' }, numeric(1))
#' summary(pseudo_z)
#' @md
#' @export
generate_pseudo_dyads <- function(
  dyad_list,
  n_pairs = NULL,
  keep_roles = TRUE
) {
  dyads <- .pseudo_extract_dyads(dyad_list)
  dyad_names <- .pseudo_dyad_names(dyad_list)
  if (!is.null(n_pairs) && !.is_count(n_pairs)) {
    cli::cli_abort("{.arg n_pairs} must be NULL or a single positive integer.")
  }
  .check_keep_roles(keep_roles)

  n <- length(dyads)
  if (keep_roles) {
    # x of dyad i with y of dyad j, i != j (i outer, j inner).
    g <- expand.grid(y_dyad = seq_len(n), x_dyad = seq_len(n))
    g <- g[g$x_dyad != g$y_dyad, ]
    pairs <- data.frame(
      x_dyad = g$x_dyad,
      x_role = "x",
      y_dyad = g$y_dyad,
      y_role = "y",
      stringsAsFactors = FALSE
    )
  } else {
    # Unordered pairs of the 2N series from different dyads, pooled as
    # rMEA's shuffle() does (1x, ..., Nx, 1y, ..., Ny); the earlier series
    # in that order becomes x, so a mixed-role pair keeps an x series as x.
    lab_dyad <- rep(seq_len(n), times = 2)
    lab_role <- rep(c("x", "y"), each = n)
    ab <- utils::combn(2 * n, 2)
    keep <- lab_dyad[ab[1, ]] != lab_dyad[ab[2, ]]
    ab <- ab[, keep, drop = FALSE]
    pairs <- data.frame(
      x_dyad = lab_dyad[ab[1, ]],
      x_role = lab_role[ab[1, ]],
      y_dyad = lab_dyad[ab[2, ]],
      y_role = lab_role[ab[2, ]],
      stringsAsFactors = FALSE
    )
  }

  n_total <- nrow(pairs)
  if (!is.null(n_pairs)) {
    if (n_pairs > n_total) {
      cli::cli_abort(
        "{.arg n_pairs} ({n_pairs}) exceeds the {n_total} possible \\
        pseudo-dyads."
      )
    }
    pairs <- pairs[sample.int(n_total, n_pairs), , drop = FALSE]
  }
  rownames(pairs) <- NULL

  n_cropped <- 0L
  out <- lapply(seq_len(nrow(pairs)), function(k) {
    sx <- dyads[[pairs$x_dyad[k]]][[pairs$x_role[k]]]
    sy <- dyads[[pairs$y_dyad[k]]][[pairs$y_role[k]]]
    m <- min(length(sx), length(sy))
    if (length(sx) != length(sy)) {
      n_cropped <<- n_cropped + 1L
    }
    list(x = sx[seq_len(m)], y = sy[seq_len(m)])
  })
  if (n_cropped > 0) {
    cli::cli_inform(
      "Cropped {n_cropped} pseudo-dyads to their shorter series \\
      (start-aligned)."
    )
  }
  name_of <- function(d) {
    if (is.null(dyad_names)) rep(NA_character_, length(d)) else dyad_names[d]
  }
  pairs <- data.frame(
    x_dyad = pairs$x_dyad,
    x_role = pairs$x_role,
    x_name = name_of(pairs$x_dyad),
    y_dyad = pairs$y_dyad,
    y_role = pairs$y_role,
    y_name = name_of(pairs$y_dyad),
    stringsAsFactors = FALSE
  )
  if (!is.null(dyad_names)) {
    names(out) <- paste0(
      pairs$x_name,
      ":",
      pairs$x_role,
      "|",
      pairs$y_name,
      ":",
      pairs$y_role
    )
  }
  attr(out, "sources") <- pairs

  out
}

# --- Pseudo-dyad helpers (internal) --------------------------------------

# Extract and type-check every dyad in a dyad_list (>= 2 dyads).
.pseudo_extract_dyads <- function(dyad_list) {
  if (is.data.frame(dyad_list)) {
    cli::cli_abort(
      "{.arg dyad_list} must be a list of dyads, not a single data frame. \\
      Wrap one data frame per dyad in a list."
    )
  }
  if (!is.list(dyad_list) || length(dyad_list) < 2) {
    cli::cli_abort(
      "{.arg dyad_list} must contain at least two dyads (a list with one \\
      element per dyad)."
    )
  }
  lapply(seq_along(dyad_list), function(i) {
    xy <- .extract_xy(dyad_list[[i]], i)
    if (!is.numeric(xy$x) || !is.numeric(xy$y)) {
      cli::cli_abort(
        "Both series of dyad {i} in {.arg dyad_list} must be numeric."
      )
    }
    if (length(xy$x) == 0 || length(xy$y) == 0) {
      cli::cli_abort(
        "Both series of dyad {i} in {.arg dyad_list} must have at least one \\
        sample."
      )
    }
    list(x = as.double(xy$x), y = as.double(xy$y))
  })
}

# The dyad names of a dyad_list: a character vector when every name is
# non-empty, not NA, and unique; NULL when the list has no names. Any other
# naming aborts, because the names are used as dyad IDs in the output.
.pseudo_dyad_names <- function(dyad_list) {
  nms <- names(dyad_list)
  if (is.null(nms)) {
    return(NULL)
  }
  if (anyNA(nms) || any(nms == "")) {
    cli::cli_abort(
      "The names of {.arg dyad_list} must all be non-empty and not \\
      {.val {NA}}, or the list must have no names."
    )
  }
  dups <- unique(nms[duplicated(nms)])
  if (length(dups) > 0) {
    cli::cli_abort(c(
      "The names of {.arg dyad_list} must be unique.",
      "x" = "{.val {dups}} {?is/are} repeated."
    ))
  }
  nms
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

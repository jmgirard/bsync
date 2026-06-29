#' Diagnose Missing Data Gaps in a Time Series
#'
#' @param x A numeric vector.
#' @return A data frame with summary statistics about missing values.
#' @export
diagnose_ts_gaps <- function(x) {
  if (!is.numeric(x)) {
    cli::cli_abort("{.arg x} must be a numeric vector.")
  }

  total_na <- sum(is.na(x))
  percent_na <- (total_na / length(x)) * 100

  na_runs <- rle(is.na(x))
  gap_lengths <- na_runs$lengths[na_runs$values == TRUE]

  max_gap <- if (length(gap_lengths) > 0) max(gap_lengths) else 0
  total_gaps <- length(gap_lengths)

  data.frame(
    total_obs = length(x),
    total_na = total_na,
    percent_na = round(percent_na, 2),
    total_gaps = total_gaps,
    max_gap_length = max_gap
  )
}

#' Impute Missing Values in Continuous Time Series with Metadata
#'
#' @param x A numeric vector representing the time series.
#' @param method Character string specifying the interpolation method ("linear" or "spline").
#' @param maxgap Integer specifying the maximum number of consecutive NAs to impute.
#'   Gaps larger than this will be left as NA.
#' @return A numeric vector with small gaps imputed, containing an "imputation_summary" attribute.
#' @export
impute_ts_gaps <- function(x, method = c("linear", "spline"), maxgap = 5) {
  method <- match.arg(method)

  if (!is.numeric(x)) {
    cli::cli_abort("{.arg x} must be a numeric vector.")
  }

  original_na <- sum(is.na(x))

  # If no missing data, attach summary and return early
  if (original_na == 0) {
    attr(x, "imputation_summary") <- list(
      method = method,
      maxgap_used = maxgap,
      values_imputed = 0,
      values_left_na = 0
    )
    return(x)
  }

  # Identify valid data points
  valid_idx <- which(!is.na(x))
  if (length(valid_idx) < 2) {
    cli::cli_warn("Not enough valid data points to perform interpolation.")
    attr(x, "imputation_summary") <- list(
      method = method,
      maxgap_used = maxgap,
      values_imputed = 0,
      values_left_na = original_na
    )
    return(x)
  }

  # Perform global imputation
  x_imp <- x

  if (method == "linear") {
    # rule = 1 ensures approx does not extrapolate
    imp_res <- stats::approx(
      x = valid_idx,
      y = x[valid_idx],
      xout = seq_along(x),
      rule = 1
    )
    x_imp <- imp_res$y
  } else if (method == "spline") {
    # Spline extrapolates by default, requiring manual cleanup below
    imp_res <- stats::spline(
      x = valid_idx,
      y = x[valid_idx],
      xout = seq_along(x),
      method = "fmm"
    )
    x_imp <- imp_res$y
  }

  # Re-apply NAs to leading and trailing edges to strictly prevent extrapolation
  first_valid <- min(valid_idx)
  last_valid <- max(valid_idx)

  if (first_valid > 1) {
    x_imp[1:(first_valid - 1)] <- NA
  }
  if (last_valid < length(x)) {
    x_imp[(last_valid + 1):length(x)] <- NA
  }

  # Re-apply NAs to internal gaps that exceed maxgap
  na_runs <- rle(is.na(x))
  long_gaps <- which(na_runs$values == TRUE & na_runs$lengths > maxgap)

  if (length(long_gaps) > 0) {
    end_idx <- cumsum(na_runs$lengths)
    start_idx <- end_idx - na_runs$lengths + 1

    for (i in long_gaps) {
      x_imp[start_idx[i]:end_idx[i]] <- NA
    }

    cli::cli_warn(
      "Found {length(long_gaps)} gap{?s} exceeding {.arg maxgap} ({maxgap}); left as NA."
    )
  }

  # Calculate post-imputation metadata
  remaining_na <- sum(is.na(x_imp))
  imputed_count <- original_na - remaining_na

  # Attach metadata as an attribute
  attr(x_imp, "imputation_summary") <- list(
    method = method,
    maxgap_used = maxgap,
    values_imputed = imputed_count,
    values_left_na = remaining_na
  )

  return(x_imp)
}

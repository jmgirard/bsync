# Main Functions ----------------------------------------------------------

#' Find Optimum (Peak or Valley) in Windowed Analyses
#'
#' @param obj An object of class "wcc_res" or "wdtw_res".
#' @param L_size An odd integer specifying the size of the local search region.
#'   Ignored if `search_method = "global"`. Default is `NULL`.
#' @param strict_monotonic Logical indicating whether to strictly enforce
#'   monotonic flanks around the extremum. Ignored if `search_method = "global"`. Default is FALSE.
#' @param find_min Logical indicating whether to search for local minima instead
#'   of local maxima. If `NULL` (the default), the function automatically
#'   searches for maxima (`FALSE`) for cross-correlation ("wcc_res") and
#'   minima (`TRUE`) for distance metrics ("wdtw_res").
#' @param search_method Character string specifying "local" or "global" search.
#'   "local" searches symmetrically outward from lag 0. "global" searches the
#'   entire window for the absolute extremum. If `NULL`, defaults to "local"
#'   for "wcc_res" and "global" for "wdtw_res".
#' @param threshold A numeric value. For WCC (`find_min = FALSE`), optima with an
#'   absolute value below this threshold are set to NA. For WDTW (`find_min = TRUE`),
#'   optima with a distance above this threshold are set to NA. Default is `NULL`.
#' @return A data frame of class "wcc_optima" or "wdtw_optima".
#' @examples
#' wcc_res <- wcc(sim_dyad$x_A, sim_dyad$x_B, window_size = 96, lag_max = 10)
#' optima <- pick_optima(wcc_res, L_size = 9)
#' head(optima)
#' @export
pick_optima <- function(
  obj,
  L_size = NULL,
  strict_monotonic = FALSE,
  find_min = NULL,
  search_method = NULL,
  threshold = NULL
) {
  if (inherits(obj, "wcc_res")) {
    metric_col <- "wcc"
    out_class <- "wcc_optima"
    if (is.null(find_min)) {
      find_min <- FALSE
    }
    if (is.null(search_method)) search_method <- "local"
  } else if (inherits(obj, "wdtw_res")) {
    metric_col <- "dtw_dist"
    out_class <- "wdtw_optima"
    if (is.null(find_min)) {
      find_min <- TRUE
    }
    if (is.null(search_method)) search_method <- "global"
  } else {
    cli::cli_abort(
      "Input {.arg obj} must be a {.cls wcc_res} or {.cls wdtw_res} object."
    )
  }

  tau_max <- obj$settings$lag_max
  df <- obj$results_df
  df <- df[order(df$i, df$tau), ]

  if (search_method == "local") {
    if (is.null(L_size)) {
      cli::cli_abort(
        "{.arg L_size} must be provided when using {.code search_method = 'local'}."
      )
    }

    i_vals <- unique(df$i)
    i_factor <- factor(df$i, levels = i_vals)
    metric_list <- split(df[[metric_col]], i_factor)

    out_df <- pick_optima_cpp(
      metric_list = metric_list,
      i_vals = i_vals,
      tau_max = tau_max,
      L_size = L_size,
      strict_monotonic = strict_monotonic,
      find_min = find_min
    )
  } else if (search_method == "global") {
    out_df <- df |>
      dplyr::filter(!is.na(.data[[metric_col]]))

    if (find_min) {
      out_df <- out_df |>
        dplyr::slice_min(
          order_by = .data[[metric_col]],
          n = 1,
          with_ties = FALSE,
          by = i
        )
    } else {
      out_df <- out_df |>
        dplyr::slice_max(
          order_by = .data[[metric_col]],
          n = 1,
          with_ties = FALSE,
          by = i
        )
    }

    out_df <- out_df |>
      dplyr::select(
        i,
        optimum_lag = tau,
        optimum_value = dplyr::all_of(metric_col)
      ) |>
      as.data.frame()
  } else {
    cli::cli_abort(
      "{.arg search_method} must be either {.val local} or {.val global}."
    )
  }

  if (!is.null(threshold)) {
    if (find_min) {
      bad_idx <- which(out_df$optimum_value > threshold)
    } else {
      bad_idx <- which(abs(out_df$optimum_value) < threshold)
    }

    out_df$optimum_lag[bad_idx] <- NA
    out_df$optimum_value[bad_idx] <- NA
  }

  attr(out_df, "search_method") <- search_method
  attr(out_df, "find_min") <- find_min
  attr(out_df, "threshold") <- threshold
  attr(out_df, "has_time") <- isTRUE(obj$settings$has_time)

  if (search_method == "local") {
    attr(out_df, "L_size") <- L_size
    attr(out_df, "strict_monotonic") <- strict_monotonic
  }

  structure(out_df, class = c(out_class, "data.frame"))
}

# S3 Methods --------------------------------------------------------------

#' Print method for wcc_optima objects
#'
#' @param x An object of class "wcc_optima".
#' @param n An integer specifying how many rows to print. Default is 5.
#' @param ... Additional arguments (not used).
#' @return Returns `x` invisibly.
#' @export
print.wcc_optima <- function(x, n = 5, ...) {
  print_optima(x, n, title = "WCC Optima Results")
}

#' Print method for wdtw_optima objects
#'
#' @param x An object of class "wdtw_optima".
#' @param n An integer specifying how many rows to print. Default is 5.
#' @param ... Additional arguments (not used).
#' @return Returns `x` invisibly.
#' @export
print.wdtw_optima <- function(x, n = 5, ...) {
  print_optima(x, n, title = "WDTW Optima Results")
}

#' Summary method for wcc_optima objects
#'
#' @param object An object of class "wcc_optima".
#' @param ... Additional arguments (not used).
#' @return Returns `object` invisibly.
#' @export
summary.wcc_optima <- function(object, ...) {
  summary_optima(object, title = "WCC Optima Summary")
}

#' Summary method for wdtw_optima objects
#'
#' @param object An object of class "wdtw_optima".
#' @param ... Additional arguments (not used).
#' @return Returns `object` invisibly.
#' @export
summary.wdtw_optima <- function(object, ...) {
  summary_optima(object, title = "WDTW Optima Summary")
}

# Internal Helpers --------------------------------------------------------

#' @noRd
print_optima <- function(x, n, title) {
  search_method <- attr(x, "search_method")
  find_min <- attr(x, "find_min")
  thresh <- attr(x, "threshold")

  total_optima <- nrow(x)
  valid_optima <- sum(!is.na(x$optimum_lag))

  cli::cli_h1(title)

  cli::cli_dl(c(
    "Total Windows Analyzed" = "{total_optima}",
    "Valid Optima Found" = "{valid_optima} ({round(valid_optima / total_optima * 100, 1)}%)",
    "Search Method" = "{search_method}",
    "Search Mode" = "{ifelse(find_min, 'Valleys (Minima)', 'Peaks (Maxima)')}",
    "Threshold Applied" = "{ifelse(is.null(thresh), 'None', thresh)}"
  ))

  if (search_method == "local") {
    cli::cli_dl(c(
      "Local Search Size" = "{attr(x, 'L_size')}",
      "Strict Monotonic" = "{attr(x, 'strict_monotonic')}"
    ))
  }

  if (total_optima > 0) {
    cli::cli_text("Showing the first {min(n, total_optima)} result{?s}:")
    print(utils::head(as.data.frame(x), n), row.names = FALSE)

    if (total_optima > n) {
      remaining <- total_optima - n
      cli::cli_text(
        "{cli::col_grey('# ... with ', remaining, ' more row', ifelse(remaining == 1, '', 's'))}"
      )
    }
  } else {
    cli::cli_alert_info("No optima found matching the criteria.")
  }

  invisible(x)
}

#' @noRd
summary_optima <- function(object, title) {
  total_optima <- nrow(object)
  valid_idx <- !is.na(object$optimum_lag)
  valid_optima <- sum(valid_idx)
  missing_optima <- total_optima - valid_optima

  cli::cli_h1(title)

  cli::cli_h2("Completeness")
  cli::cli_bullets(c(
    "*" = "Total time windows: {total_optima}",
    "*" = "Valid optima retained: {valid_optima} ({round(valid_optima / total_optima * 100, 1)}%)",
    "*" = "Optima dropped (NA): {missing_optima} ({round(missing_optima / total_optima * 100, 1)}%)"
  ))

  if (valid_optima > 0) {
    lags <- object$optimum_lag[valid_idx]
    vals <- object$optimum_value[valid_idx]

    pos_lags <- sum(lags > 0)
    neg_lags <- sum(lags < 0)
    zero_lags <- sum(lags == 0)

    cli::cli_h2("Lag Directionality (Leadership)")
    cli::cli_bullets(c(
      "*" = "Positive Lags (x leads y): {pos_lags} ({round(pos_lags / valid_optima * 100, 1)}%)",
      "*" = "Negative Lags (y leads x): {neg_lags} ({round(neg_lags / valid_optima * 100, 1)}%)",
      "*" = "Zero Lags (Simultaneous):  {zero_lags} ({round(zero_lags / valid_optima * 100, 1)}%)"
    ))

    cli::cli_h2("Optimum Value Distribution")
    q_vals <- stats::quantile(vals, probs = c(0, 0.25, 0.5, 0.75, 1))
    print(round(q_vals, 4))
  }

  invisible(object)
}

# Tidy interface for bsync_surface objects (M5) --------------------------------
#
# Implements generics::tidy(), generics::glance(), and tibble::as_tibble() for
# the "bsync_surface" superclass so users can use the standard tidyverse verbs.
#
# Shape conventions:
#   tidy()     — one row per cell of results_df (the full windowed surface)
#   glance()   — one row per estimator run (aggregate scalar(s) + settings)
#   as_tibble()— same as tidy() but via the tibble generic
#
# Granger note: two aggregate statistics (f_xy, f_yx) are stacked as separate
# columns in glance(), and as a pair of rows (by direction) in an auxiliary
# column rather than pivoting the surface (which already has f_xy / f_yx cols).

# Register generics::tidy and generics::glance
# (generics is in Imports, so @importFrom is needed for NAMESPACE)

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom tibble as_tibble
#' @export
tibble::as_tibble


# tidy() ----------------------------------------------------------------------

#' Tidy a bsync_surface object into a tibble of per-cell results
#'
#' Returns one row per cell in `results_df`: window position × lag (or just
#' window position for Granger). Column names match the underlying estimator.
#'
#' @param x A `bsync_surface` object (`wcc_res`, `wdtw_res`, or `wgranger_res`).
#' @param ... Additional arguments (not used).
#' @return A [tibble::tibble()].
#' @seealso [glance.bsync_surface()], [as_tibble.bsync_surface()]
#' @exportS3Method generics::tidy bsync_surface
tidy.bsync_surface <- function(x, ...) {
  tibble::as_tibble(x$results_df)
}


# glance() --------------------------------------------------------------------

#' One-row summary of a bsync_surface object
#'
#' Returns a single-row tibble with the aggregate statistic(s) and key settings
#' so multiple estimator runs can be compared with [dplyr::bind_rows()].
#'
#' For WCC/WDTW the aggregate is a single named column (`mean_abs_z`, `peak`,
#' or `mean_distance`). For Granger, `f_xy` and `f_yx` are separate columns.
#'
#' @param x A `bsync_surface` object (`wcc_res`, `wdtw_res`, or `wgranger_res`).
#' @param ... Additional arguments (not used).
#' @return A one-row [tibble::tibble()].
#' @seealso [tidy.bsync_surface()], [as_tibble.bsync_surface()]
#' @exportS3Method generics::glance bsync_surface
glance.bsync_surface <- function(x, ...) {
  agg <- as.list(x$aggregate) # names preserved; scalar or two-element list
  cfg <- x$settings

  # n_windows is not stored in settings; derive from the surface
  n_unique_i <- if (!is.null(x$results_df$i)) {
    length(unique(x$results_df$i))
  } else {
    NA_integer_
  }

  # Common fields available on every surface
  common <- list(
    n_windows        = n_unique_i,
    window_size      = cfg$window_size %||% NA_integer_,   # rlang::`%||%`
    window_increment = cfg$window_increment %||% NA_integer_
  )

  # Estimator-specific settings from $settings (keys vary by estimator)
  extra <- list()
  if (!is.null(cfg$lag_max)) extra$lag_max <- cfg$lag_max
  if (!is.null(cfg$lag_increment)) extra$lag_increment <- cfg$lag_increment
  if (!is.null(cfg$statistic)) extra$statistic <- cfg$statistic
  if (!is.null(cfg$scale_method)) extra$scale_method <- cfg$scale_method
  if (!is.null(cfg$distance_metric)) extra$distance_metric <- cfg$distance_metric
  if (!is.null(cfg$ar_order)) extra$ar_order <- cfg$ar_order

  tibble::as_tibble(c(agg, common, extra))
}


# as_tibble() -----------------------------------------------------------------

#' Convert a bsync_surface object to a tibble
#'
#' Alias for [tidy.bsync_surface()]: returns one row per cell of `results_df`.
#'
#' @param x A `bsync_surface` object.
#' @param ... Additional arguments (not used).
#' @return A [tibble::tibble()].
#' @seealso [tidy.bsync_surface()], [glance.bsync_surface()]
#' @exportS3Method tibble::as_tibble bsync_surface
as_tibble.bsync_surface <- function(x, ...) {
  tibble::as_tibble(x$results_df)
}

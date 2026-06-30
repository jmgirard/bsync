#' Calculate 1D Velocity
#'
#' Calculates the velocity (rate of change including direction) along a single axis using finite difference methods.
#'
#' @param t A numeric vector representing time.
#' @param x A numeric vector of coordinates.
#' @param n An integer specifying the step size for the difference calculations. Default is 1.
#' @param method A character string specifying the method to use: "central", "forward", or "backward". Default is "central".
#' @param fill_edges A logical indicating whether to automatically use forward and backward differences to estimate velocities at the boundaries when `method = "central"`. Default is `TRUE`.
#' @return A numeric vector of velocities the same length as the input vectors.
#' @examples
#' # Signed 1D velocity of one coordinate over time
#' v <- calc_velocity_1d(t = sim_dyad$time, x = sim_dyad$x_A)
#' head(v)
#' @export
calc_velocity_1d <- function(
  t,
  x,
  n = 1,
  method = c("central", "forward", "backward"),
  fill_edges = TRUE
) {
  method <- match.arg(method)

  calc_fwd <- function() (dplyr::lead(x, n) - x) / (dplyr::lead(t, n) - t)
  calc_bwd <- function() (x - dplyr::lag(x, n)) / (t - dplyr::lag(t, n))

  if (method == "forward") {
    return(calc_fwd())
  }
  if (method == "backward") {
    return(calc_bwd())
  }

  # Default to central difference
  dt <- dplyr::lead(t, n) - dplyr::lag(t, n)
  dx <- dplyr::lead(x, n) - dplyr::lag(x, n)
  v <- dx / dt

  if (fill_edges) {
    v <- dplyr::coalesce(v, calc_fwd(), calc_bwd())
  }

  v
}

#' Calculate 1D Speed
#'
#' Calculates the absolute speed along a single axis (e.g., x-axis) using finite difference methods.
#'
#' @param t A numeric vector representing time.
#' @param x A numeric vector of coordinates.
#' @param n An integer specifying the step size for the difference calculations. Default is 1.
#' @param method A character string specifying the method to use: "central", "forward", or "backward". Default is "central".
#' @param fill_edges A logical indicating whether to automatically use forward and backward differences to estimate speeds at the boundaries when `method = "central"`. Default is `TRUE`.
#' @return A numeric vector of speeds the same length as the input vectors.
#' @examples
#' # Unsigned 1D speed (absolute rate of change)
#' s <- calc_speed_1d(t = sim_dyad$time, x = sim_dyad$x_A)
#' head(s)
#' @export
calc_speed_1d <- function(
  t,
  x,
  n = 1,
  method = c("central", "forward", "backward"),
  fill_edges = TRUE
) {
  method <- match.arg(method)

  calc_fwd <- function() abs((dplyr::lead(x, n) - x) / (dplyr::lead(t, n) - t))
  calc_bwd <- function() abs((x - dplyr::lag(x, n)) / (t - dplyr::lag(t, n)))

  if (method == "forward") {
    return(calc_fwd())
  }
  if (method == "backward") {
    return(calc_bwd())
  }

  # Default to central difference
  dt <- dplyr::lead(t, n) - dplyr::lag(t, n)
  dx <- dplyr::lead(x, n) - dplyr::lag(x, n)
  v <- abs(dx / dt)

  if (fill_edges) {
    v <- dplyr::coalesce(v, calc_fwd(), calc_bwd())
  }

  v
}

#' Calculate 2D Speed
#'
#' Calculates the magnitude of the 2D velocity vector (speed) using finite difference methods.
#'
#' @param t A numeric vector representing time.
#' @param x A numeric vector of x-coordinates.
#' @param y A numeric vector of y-coordinates.
#' @param n An integer specifying the step size for the difference calculations. Default is 1.
#' @param method A character string specifying the method to use: "central", "forward", or "backward". Default is "central".
#' @param fill_edges A logical indicating whether to automatically use forward and backward differences to estimate speeds at the boundaries when `method = "central"`. Default is `TRUE`.
#' @return A numeric vector of speeds the same length as the input vectors.
#' @examples
#' # 2D speed from a pair of coordinate channels
#' s <- calc_speed_2d(t = sim_dyad$time, x = sim_dyad$x_A, y = sim_dyad$y_A)
#' head(s)
#' @export
calc_speed_2d <- function(
  t,
  x,
  y,
  n = 1,
  method = c("central", "forward", "backward"),
  fill_edges = TRUE
) {
  method <- match.arg(method)

  calc_fwd <- function() {
    dt <- dplyr::lead(t, n) - t
    dx <- dplyr::lead(x, n) - x
    dy <- dplyr::lead(y, n) - y
    sqrt((dx / dt)^2 + (dy / dt)^2)
  }

  calc_bwd <- function() {
    dt <- t - dplyr::lag(t, n)
    dx <- x - dplyr::lag(x, n)
    dy <- y - dplyr::lag(y, n)
    sqrt((dx / dt)^2 + (dy / dt)^2)
  }

  if (method == "forward") {
    return(calc_fwd())
  }
  if (method == "backward") {
    return(calc_bwd())
  }

  # Default to central difference
  dt <- dplyr::lead(t, n) - dplyr::lag(t, n)
  dx <- dplyr::lead(x, n) - dplyr::lag(x, n)
  dy <- dplyr::lead(y, n) - dplyr::lag(y, n)
  v <- sqrt((dx / dt)^2 + (dy / dt)^2)

  if (fill_edges) {
    v <- dplyr::coalesce(v, calc_fwd(), calc_bwd())
  }

  v
}

#' Calculate 3D Speed
#'
#' Calculates the magnitude of the 3D velocity vector (speed) using finite difference methods.
#'
#' @param t A numeric vector representing time.
#' @param x A numeric vector of x-coordinates.
#' @param y A numeric vector of y-coordinates.
#' @param z A numeric vector of z-coordinates.
#' @param n An integer specifying the step size for the difference calculations. Default is 1.
#' @param method A character string specifying the method to use: "central", "forward", or "backward". Default is "central".
#' @param fill_edges A logical indicating whether to automatically use forward and backward differences to estimate speeds at the boundaries when `method = "central"`. Default is `TRUE`.
#' @return A numeric vector of speeds the same length as the input vectors.
#' @examples
#' # 3D speed from x/y/z coordinate channels
#' s <- calc_speed_3d(
#'   t = sim_dyad$time, x = sim_dyad$x_A, y = sim_dyad$y_A, z = sim_dyad$z_A
#' )
#' head(s)
#' @export
calc_speed_3d <- function(
  t,
  x,
  y,
  z,
  n = 1,
  method = c("central", "forward", "backward"),
  fill_edges = TRUE
) {
  method <- match.arg(method)

  calc_fwd <- function() {
    dt <- dplyr::lead(t, n) - t
    dx <- dplyr::lead(x, n) - x
    dy <- dplyr::lead(y, n) - y
    dz <- dplyr::lead(z, n) - z
    sqrt((dx / dt)^2 + (dy / dt)^2 + (dz / dt)^2)
  }

  calc_bwd <- function() {
    dt <- t - dplyr::lag(t, n)
    dx <- x - dplyr::lag(x, n)
    dy <- y - dplyr::lag(y, n)
    dz <- z - dplyr::lag(z, n)
    sqrt((dx / dt)^2 + (dy / dt)^2 + (dz / dt)^2)
  }

  if (method == "forward") {
    return(calc_fwd())
  }
  if (method == "backward") {
    return(calc_bwd())
  }

  # Default to central difference
  dt <- dplyr::lead(t, n) - dplyr::lag(t, n)
  dx <- dplyr::lead(x, n) - dplyr::lag(x, n)
  dy <- dplyr::lead(y, n) - dplyr::lag(y, n)
  dz <- dplyr::lead(z, n) - dplyr::lag(z, n)
  v <- sqrt((dx / dt)^2 + (dy / dt)^2 + (dz / dt)^2)

  if (fill_edges) {
    v <- dplyr::coalesce(v, calc_fwd(), calc_bwd())
  }

  v
}

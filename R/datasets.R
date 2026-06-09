#' Generate a reverse-L necessity scatter.
#'
#' The condition X is uniform on [0, 100]; the realised outcome is a fraction
#' \code{U^k} of the linear ceiling \code{ci + cs * X}, plus Gaussian noise. This
#' is the canonical scatter NCA and QNCA detect: a full lower-right triangle and
#' an empty upper-left corner.
#'
#' @param n number of cases.
#' @param ci,cs intercept and slope of the deterministic ceiling.
#' @param k efficiency exponent (larger pushes cases toward the floor).
#' @param noise standard deviation of additive Gaussian noise.
#' @return data frame with columns X and Y, both clamped to [0, 100].
#' @export
generate_reverse_L <- function(n, ci = 15, cs = 0.85, k = 2, noise = 4) {
  X <- stats::runif(n, 0, 100)
  ceiling <- ci + cs * X
  eff <- stats::runif(n)^k
  Y <- eff * ceiling + stats::rnorm(n, 0, noise)
  data.frame(X = pmin(pmax(X, 0), 100), Y = pmin(pmax(Y, 0), 100))
}

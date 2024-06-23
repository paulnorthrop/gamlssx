#' GEV Distribution Expected Information
#'
#' Calculates the expected information matrix for the GEV distribution.
#'
#' @param scale,shape Numeric vectors. Respective values of the GEV parameters
#'   scale parameter \eqn{\sigma} and shape parameter \eqn{\xi}. For
#'   `gevExpInfo`, `scale` and `shape` must have length 1.
#' @param eps A numeric scalar. For values of \eqn{\xi} in `shape` that lie in
#'   `(-eps, eps)` an approximation is used instead of a direct calculation.
#'   See **Details**. If `eps` is a vector then only the first element is used.
#' @details `gevExpInfo` calculates, for single pair of values
#'   \eqn{(\sigma, \xi) = } `(scale, shape)`, the expected information matrix for a
#'   single observation from a GEV distribution with distribution function
#'   \deqn{F(x) = P(X \leq x) = \exp\left\{ -\left[ 1+\xi\left(\frac{x-\mu}{\sigma}\right)
#'   \right]_+^{-1/\xi} \right\},}
#'   where \eqn{x_+ = \max(x, 0)}.
#'   The GEV expected information is defined only for \eqn{\xi > -0.5} and does
#'   not depend on the value of \eqn{\mu}.
#'
#'   The other functions are vectorized and calculate the individual
#'   contributions to the expected information matrix. For example, `gev11e`
#'   calculates the expectation \eqn{i_{\mu\mu}} of the negated second
#'   derivative of the GEV log-density with respect to \eqn{\mu}, that is, each
#'   `1` indicates one derivative with respect to \eqn{\mu}. Similarly, `2`
#'   denotes one derivative with respect to \eqn{\sigma} and `3` one derivative
#'   with respect to \eqn{\xi}, so that, for example, `gev23e` calculates the
#'   expectation \eqn{i_{\sigma\xi}} of the negated GEV log-density after one
#'   taking one derivative with respect to \eqn{\sigma} and one derivative with
#'   respect to \eqn{\xi}. Note that \eqn{i_{\xi\xi}}, calculated using
#'   `gev33e`, depends only on \eqn{\xi}.
#'
#'   The expectation in `gev11e` can be calculated in a direct way for all
#'   \eqn{\xi > -0.5}. For the other components, direct calculation of the
#'   expectation is unstable when \eqn{\xi} is close to 0. Instead, we use
#'   a quadratic approximation over `(-eps, eps)`, from a Lagrangian
#'   interpolation of the values from the direct calculation for \eqn{\xi = }
#'   `-eps`, \eqn{0} and `eps`.
#' @returns `gevExpInfo` returns a 3 by 3 numeric matrix with row and column
#'   named `loc, scale, shape`. The other functions return a numeric vector of
#'   length equal to the maximum of the lengths of the arguments, excluding
#'   `eps`.
#' @examples
#' # Expected information matrices for ...
#' # ... scale = 2 and shape = -0.4
#' gevExpInfo(2, -0.4)
#' # ... scale = 3 and shape = 0.001
#' gevExpInfo(3, 0.001)
#' # ... scale = 3 and shape = 0
#' gevExpInfo(3, 0)
#' # ... scale = 1 and shape = 0.1
#' gevExpInfo(1, 0.1)
#'
#' # The individual components of the latter matrix
#' gev11e(1, 0.1)
#' gev12e(1, 0.1)
#' gev13e(1, 0.1)
#' gev22e(1, 0.1)
#' gev23e(1, 0.1)
#' gev33e(0.1)
#' @name gevExpInfo
NULL
## NULL

#' @rdname gevExpInfo
#' @export
gev11e <- function(scale, shape) {
  m <- max(length(scale), length(shape))
  return(rep_len(pxi(shape), m) / rep_len(scale ^ 2, m))
}

#' @rdname gevExpInfo
#' @export
gev22e <- function(scale, shape, eps = 3e-3) {
  eps <- eps[1]
  val <- gevExpInfoComp(fun = gev22eFn, fun0 = gev22e0Constant, xi = shape,
                        eps = eps)
  m <- max(length(scale), length(shape))
  return(rep_len(val, m) / rep_len(scale ^ 2, m))
}

#' @rdname gevExpInfo
#' @export
gev33e <- function(shape, eps = 3e-3) {
  eps <- eps[1]
  val <- gevExpInfoComp(fun = gev33eFn, fun0 = gev33e0Constant, xi = shape,
                        eps = eps)
  return(val)
}

#' @rdname gevExpInfo
#' @export
gev12e <- function(scale, shape, eps = 3e-3) {
  eps <- eps[1]
  val <- gevExpInfoComp(fun = gev12eFn, fun0 = gev12e0Constant, xi = shape,
                        eps = eps)
  m <- max(length(scale), length(shape))
  return(rep_len(val, m) / rep_len(scale ^ 2, m))
}

#' @rdname gevExpInfo
#' @export
gev13e <- function(scale, shape, eps = 3e-3) {
  eps <- eps[1]
  val <- gevExpInfoComp(fun = gev13eFn, fun0 = gev13e0Constant, xi = shape,
                        eps = eps)
  m <- max(length(scale), length(shape))
  return(rep_len(val, m) / rep_len(scale, m))
}

#' @rdname gevExpInfo
#' @export
gev23e <- function(scale, shape, eps = 3e-3) {
  eps <- eps[1]
  val <- gevExpInfoComp(fun = gev23eFn, fun0 = gev23e0Constant, xi = shape,
                        eps = eps)
  m <- max(length(scale), length(shape))
  return(rep_len(val, m) / rep_len(scale, m))
}

#' @rdname gevExpInfo
#' @export
gevExpInfo <- function(scale, shape, eps = 3e-3) {
  if (shape <= -0.5) {
    stop("The GEV expected information is undefined for shape <= -0.5")
  }
  # The expected information does not depend on loc
  val <- matrix(NA, 3, 3)
  val[1, 1] <- gev11e(scale, shape)
  val[2, 2] <- gev22e(scale, shape, eps)
  val[3, 3] <- gev33e(shape, eps)
  val[2, 1] <- val[1, 2] <- gev12e(scale, shape, eps)
  val[3, 1] <- val[1, 3] <- gev13e(scale, shape, eps)
  val[3, 2] <- val[2, 3] <- gev23e(scale, shape, eps)
  dimnames(val) <- list(c("loc", "scale", "shape"), c("loc", "scale", "shape"))
  return(val)
}

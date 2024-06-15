#' Internal gamlssx functions
#'
#' Internal gamlssx functions
#' @details
#' These functions are not intended to be called by the user.
#' @name gamlssx-internal
#' @keywords internal
NULL

# =============== For calculating the GEV expected information ============== #

# --------------------- Constants and helper functions ---------------------- #

#' @keywords internal
#' @rdname evils-internal
EulersConstant <- 0.57721566490153286060651209008240243104215933593992

#' @keywords internal
#' @rdname evils-internal
AperysConstant <- 1.202056903159594285399738161511449990764986292

#' @keywords internal
#' @rdname evils-internal
pxi <- function(xi) {
  return((1 + xi) ^ 2 * gamma(1 + 2 * xi))
}

#' @keywords internal
#' @rdname evils-internal
qxi <- function(xi) {
  return(gamma(2 + xi) * (digamma(1 + xi) + (1 + xi) / xi))
}

# ------------------------------ (sigma, sigma) ----------------------------- #

#' @keywords internal
#' @rdname evils-internal
gev22e0Fn <- function() {
  return(pi ^ 2 / 6 + (1 - EulersConstant) ^ 2)
}

#' @keywords internal
#' @rdname evils-internal
gev22e0Constant <- gev22e0Fn()

#' @keywords internal
#' @rdname evils-internal
gev22eFn <- function(xi) {
  return((1 - 2 * gamma(2 + xi) + pxi(xi)) / (xi ^ 2))
}

# --------------------------------- (xi, xi) -------------------------------- #

#' @keywords internal
#' @rdname evils-internal
gev33e0Fn <- function() {
  val <- pi ^ 2 / 6 - pi ^ 2 * EulersConstant / 2 + EulersConstant ^ 2 -
    EulersConstant ^ 3 - 2 * AperysConstant +
    2 * EulersConstant * AperysConstant + pi ^ 2 * EulersConstant ^ 2 / 4 +
    EulersConstant ^ 4 / 4 + 3 * pi ^ 4 / 80
  return(val)
}

#' @keywords internal
#' @rdname evils-internal
gev33e0Constant <- gev33e0Fn()

#' @keywords internal
#' @rdname evils-internal
gev33eFn <- function(xi) {
  val <- (pi ^ 2 / 6 + (1 - EulersConstant + 1 / xi) ^ 2 - 2 * qxi(xi) / xi +
            pxi(xi) / xi ^ 2) / (xi ^ 2)
  return(val)
}

# -------------------------------- (mu, sigma) ------------------------------ #

#' @keywords internal
#' @rdname evils-internal
gev12e0Fn <- function() {
  return(EulersConstant - 1)
}

#' @keywords internal
#' @rdname evils-internal
gev12e0Constant <- gev12e0Fn()

#' @keywords internal
#' @rdname evils-internal
gev12eFn <- function(xi) {
  val <- (gamma(2 + xi) - pxi(xi)) / xi
  return(val)
}

# --------------------------------- (mu, xi) -------------------------------- #

#' @keywords internal
#' @rdname evils-internal
gev13e0Fn <- function() {
  return((pi ^ 2 / 6 + EulersConstant ^ 2 - 2 * EulersConstant) / 2)
}

#' @keywords internal
#' @rdname evils-internal
gev13e0Constant <- gev13e0Fn()

#' @keywords internal
#' @rdname evils-internal
gev13eFn <- function(xi) {
  val <- (pxi(xi) / xi - qxi(xi)) / xi
  return(val)
}

# ------------------------------- (sigma, xi) ------------------------------- #

#' @keywords internal
#' @rdname evils-internal
gev23e0Fn <- function() {
  val <- (4 * EulersConstant + 4 * AperysConstant + pi ^ 2 * EulersConstant +
     2 * EulersConstant ^ 3 - pi ^ 2 - 6 * EulersConstant ^ 2) / 4
  return(val)
}

#' @keywords internal
#' @rdname evils-internal
gev23e0Constant <- gev23e0Fn()

#' @keywords internal
#' @rdname evils-internal
gev23eFn <- function(xi) {
  val <- -(1 - EulersConstant + (1 - gamma(2 + xi)) / xi - qxi(xi) +
             pxi(xi) / xi) / (xi ^ 2)
  return(val)
}

# -- Calculate a component using a quadratic approximation if xi is near 0 -- #

#' @keywords internal
#' @rdname evils-internal
gevExpInfoComp <- function(fun, fun0, xi, eps = 3e-3) {
  eps <- abs(eps)
  val <- xi
  if (any(xiNearZero <- abs(xi) < eps)) {
    aa <- fun0
    yp <- fun(eps)
    ym <- fun(-eps)
    ff <- lagrangianInterpolation(c(-eps, 0, eps), c(ym, aa, yp))
    val[xiNearZero] <- ff(xi[xiNearZero])
  }
  val[!xiNearZero] <- fun(xi[!xiNearZero])
  return(val)
}

#' @keywords internal
#' @rdname evils-internal
lagrangianInterpolation <- function(x0, y0) {
  f <- function(x) {
    sum(y0 * sapply(seq_along(x0), \(j) {
      prod(x - x0[-j])/prod(x0[j] - x0[-j])
    }))
  }
  return(Vectorize(f, "x"))
}

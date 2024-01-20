#' Generalized extreme value family distribution for fitting a GAMLSS
#'
#' The function `GEV` defines the generalized extreme value family
#' distribution, a three parameter distribution, for a
#' [`gamlss.family`][`gamlss.dist::gamlss.family`] object to be used in GAMLSS
#' fitting using the function [`gamlss`][`gamlss::gamlss`]. The functions
#' `dGEV`, `pGEV`, `qGEV` and `rGEV` define the density, distribution function,
#' quantile function and random generation for the specific parameterization of
#' the generalized extreme value distribution given in details below.
#'
#' @param mu.link Defines the `mu.link`, with `"identity"` link as the default
#' for the `mu` parameter.
#' @param sigma.link Defines the `sigma.link`, with `"log"` link as the default
#' for the `sigma` parameter.
#' @param nu.link Defines the `nu.link`, with `"identity"` link as the default
#' for the `nu` parameter.
#' @param x,q Vector of quantiles.
#' @param mu,sigma,nu Vectors of location, scale and shape parameter values.
#' @param log,log.p Logical. If `TRUE`, probabilities `eqn{p}` are given as
#'   \eqn{\log(p)}.
#' @param lower.tail Logical. If `TRUE` (the default), probabilities are
#'   \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}.
#' @param p Vector of probabilities.
#' @param n Number of observations. If `length(n) > 1`, the length is taken to
#'   be the number required.
#'
#' @details Add details.
#'
#' * Refer to Chapter 3 of Coles (2001) and Jenkinson (1955).
#' * GEV(\eqn{\mu, \sigma, \xi}). \eqn{\nu = \xi}.
#' * Give the pdf and/or cdf.
#' * Explain the initial estimates.
#' * Add code for the mean and variance.
#'
#' @return `GEV()` returns a [`gamlss.family`][`gamlss.dist::gamlss.family`]
#'   object which can be used to fit a generalized extreme value distribution
#'   in the [`gamlss`][`gamlss::gamlss`] function. `dGEV()` gives the density,
#'   `pGEV()` gives the distribution function, `qGEV()` gives the quantile
#'   function, and `rGEV()` generates random deviates.
#' @seealso [`gamlss.family`][`gamlss.dist::gamlss.family`],
#'   [`gamlss`][`gamlss::gamlss`]
#' @examples
#' # Simulate some data
#' set.seed(17012023)
#' n <- 100
#' x <- stats::runif(n)
#' mu <- 1 + 2 * x
#' sigma <- 1
#' xi <- 0.25
#' y <- nieve::rGEV(n = 1, loc = mu, scale = sigma, shape = xi)
#' plot(x, y)
#' # Fit model
#' data <- data.frame(y = as.numeric(y), x = x)
#' library(gamlss)
#' mod <- gamlss(y ~ pb(x), family = GEV, data = data)
#' plot(mod)
#' plot(data$x, data$y)
#' lines(data$x, fitted(mod))
#'
#' # Converges
#' modRS <- gamlss(y ~ pb(x), family = GEV, data = data, method = RS())
#' # Throws an error: parameters out-of-bounds leads to Inf deviance
#' modCG <- gamlss(y ~ pb(x), family = GEV, data = data, method = CG())
#' # Can we avoid the problem by halving the step lengths?
#' modCG <- gamlss(y ~ pb(x), family = GEV, data = data, method = CG(),
#'  control = gamlss.control(mu.step = 0.5, sigma.step = 0.5, nu.step = 0.5))
#  # 2 iterations of RS before switching to CG results in convergence
#' modMixed <- gamlss(y ~ pb(x), family = GEV, data = data, method = mixed())
#' @name GEV
NULL
## NULL

#' @rdname GEV
#' @export
GEV <- function(mu.link = "identity", sigma.link = "log",
                nu.link = "identity") {

  mstats <- gamlss.dist::checklink("mu.link", "GEV", substitute(mu.link),
                                   c("1/mu^2", "log", "identity"))
  dstats <- gamlss.dist::checklink("sigma.link", "GEV", substitute(sigma.link),
                                   c("inverse", "log", "identity"))
  vstats <- gamlss.dist::checklink("nu.link", "GEV",substitute(nu.link),
                                   c("inverse", "log", "identity"))

  structure(
    list(family = c("GEV", "Generalized Extreme Value"),
         parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE),
              nopar = 3,
               type = "Continuous",
            mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
            nu.link = as.character(substitute(nu.link)),
         mu.linkfun = mstats$linkfun,
      sigma.linkfun = dstats$linkfun,
         nu.linkfun = vstats$linkfun,
         mu.linkinv = mstats$linkinv,
      sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
              mu.dr = mstats$mu.eta,
           sigma.dr = dstats$mu.eta,
              nu.dr = vstats$mu.eta,
               dldm = function(y, mu, sigma, nu) {
                 dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                                   log = TRUE, deriv = TRUE)
                 dldm <- attr(dl, "gradient")[, "loc"]
                 return(dldm)
               },
             d2ldm2 = function(y, mu, sigma, nu) {
               dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                                 log = TRUE, deriv = TRUE)
               dldm <- attr(dl, "gradient")[, "loc"]
               dldm2 <- -dldm * dldm
               return(dldm2)
             },
               dldd = function(y, mu, sigma, nu) {
                 dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                                   log = TRUE, deriv = TRUE)
                 dldd <- attr(dl, "gradient")[, "scale"]
                 return(dldd)
               },
             d2ldd2 = function(y, mu, sigma, nu) {
               dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                                 log = TRUE, deriv = TRUE)
               dldd <- attr(dl, "gradient")[, "scale"]
               dldd2 <- -dldd * dldd
               return(dldd2)
             },
               dldv = function(y, mu, sigma, nu) {
                 dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                                   log = TRUE, deriv = TRUE)
                dldv <- attr(dl, "gradient")[, "shape"]
                return(dldv)
             },
             d2ldv2 = function(y, mu, sigma, nu) {
               dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                                 log = TRUE, deriv = TRUE)
               dldv <- attr(dl, "gradient")[, "shape"]
               dldv2 <- -dldv * dldv
               return(dldv2)
             },
            d2ldmdd = function(y, mu, sigma, nu) {
              dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                                log = TRUE, deriv = TRUE)
              dldm <- attr(dl, "gradient")[, "loc"]
              dldd <- attr(dl, "gradient")[, "scale"]
              dldmdd <- -dldm * dldd
              return(dldmdd)
            },
            d2ldmdv = function(y, mu, sigma, nu) {
              dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                                log = TRUE, deriv = TRUE)
              dldm <- attr(dl, "gradient")[, "loc"]
              dldv <- attr(dl, "gradient")[, "shape"]
              dldmdv <- -dldm * dldv
              return(dldmdv)
            },
            d2ldddv = function(y, mu, sigma, nu) {
              dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                          log = TRUE, deriv = TRUE)
              dldd <- attr(dl, "gradient")[, "scale"]
              dldv <- attr(dl, "gradient")[, "shape"]
              dldddv <- -dldd * dldv
              return(dldddv)
            },
#              d2ldmdd = function(y,mu,sigma,nu)  rep(0,length(y)),
#              d2ldmdv = function(y,mu,sigma,nu)  rep(0,length(y)),
#              d2ldddv = function(y,mu,sigma,nu)  rep(0,length(y)),
        G.dev.incr  = function(y, mu, sigma, nu,...) {
          val <- -2 * dGEV(x = y, mu = mu, sigma = sigma, nu = nu, log = TRUE)
#          if (any(whichInf <- is.infinite(val))) {
#            val[whichInf] <- 1e30
#          }
#          print(summary(val))
          return(val)
        },
              rqres = expression(rqres(pfun = "pGEV", type = "Continuous",
                                       y = y, mu = mu, sigma = sigma, nu = nu)),
         mu.initial = expression(mu <- y + 0.45 * sd(y)),
      sigma.initial = expression(sigma <- rep(0.78 * sd(y), length(y))),
         nu.initial = expression(nu <- rep(0.1, length(y))),
           mu.valid = function(mu) TRUE,
        sigma.valid = function(sigma) all(sigma > 0),
#           nu.valid = function(nu) TRUE,
           nu.valid = function(nu) all(mu > -0.5),
            y.valid = function(y) TRUE
    ),
  class = c("gamlss.family","family")
  )
}

#' @rdname GEV
#' @export
dGEV <- function(x, mu = 0, sigma = 1, nu = 0, log = FALSE) {
  return(nieve::dGEV(x = x, loc = mu, scale = sigma, shape = nu,
                     log = log))
}

#' @rdname GEV
#' @export
pGEV <- function(q, mu = 0, sigma = 1, nu = 0, lower.tail = TRUE,
                 log.p = FALSE) {
  return(nieve::pGEV(q = q, loc = mu, scale = sigma, shape = nu))
}

#' @rdname GEV
#' @export
qGEV <- function(p, mu = 0, sigma = 1, nu = 0, lower.tail = TRUE,
                 log.p = FALSE) {
  return(nieve::qGEV(p = p, loc = mu, scale = sigma, shape = nu))
}

#' @rdname GEV
#' @export
rGEV <- function(n, mu = 0, sigma = 1, nu = 0) {
  return(nieve::rGEV(n = n, loc = mu, scale = sigma, shape = nu))
}

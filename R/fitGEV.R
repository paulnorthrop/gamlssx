#' Fit a Generalized extreme value (GEV) GAMLSS
#'
#' Describe  [`GEV`]
#'
#' @inheritParams gamlss::gamlss
#'
#' @param scoring A character scalar. If `scoring = "fisher"` then the weights
#'   used in the fitting algorithm are based on the expected Fisher
#'   information, that is, a Fisher's scoring algorithm is used.
#'   If `scoring = "quasi` then these weights are based on the cross products
#'   of the first derivatives of the log-likelihood, leading to a quasi Newton
#'   scoring algorithm.
#' @param mu.link,sigma.link,xi.link Character scalars to set the respective
#'   link functions for the location, scale and shape parameters.
#' @param stepLength A numeric vector containing positive values. The initial
#'    values of the step lengths `mu.step`, `sigma.step` and `nu.step` passed to
#'   [`gamlss.control`][`gamlss::gamlss.control`] in the first attempt to fit
#'   the model by calling [`gamlss`][`gamlss::gamlss`]. If `stepLength` has a
#'   length that is less than 3 then `stepLength` is recycled to have length 3.
#' @param stepAttempts A non-negative integer. If the first call to
#'   [`gamlss`][`gamlss::gamlss`] throws an error then we make `stepAttempts`
#'   further attempts to fit the model, each time dividing by 2 the values
#'   of `mu.step`, `sigma.step` and `nu.step` supplied to
#'   [`gamlss.control`][`gamlss::gamlss.control`]. If `stepAttempts < 1` then
#'   no further attempts are made.
#' @param stepReduce A number greater than 1. The factor by which the step
#'   lengths in `stepLength` are reduced for each extra attempt to fit the
#'   model. The default, `stepReduce = 2` means that the step lengths are
#'   halved for each extra attempt.
#' @param eps Argument `eps` passed to [`gevExpInfo`].
#' @param ... Further arguments passed to [`gamlss`][`gamlss::gamlss`], in
#'   particular `method`, with options `RS()`, `CG()` or `mixed()`.
#'
#' @details Add details. Explain `stepAttempts` in more detail.
#'
#' @return Returns a gamlss object. See the **Value** section of
#'   [`gamlss`][`gamlss::gamlss`]. The class of the returned object is
#'   `c("gamlssx", "gamlss", "gam", "glm", "lm")`.
#' @seealso [`GEV`], [`gamlss.family`][`gamlss.dist::gamlss.family`],
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
#' data <- data.frame(y = as.numeric(y), x = x)
#' library(gamlss)
#'
#' # Fit model using the default RS method
#' mod <- fitGEV(y ~ pb(x), data = data)
#'
#' plot(mod)
#' plot(data$x, data$y)
#' lines(data$x, fitted(mod))
#' fitGEV(y ~ pb(x), data = data)
#'
#' # Fit model using the mixed method
#' mod <- fitGEV(y ~ pb(x), data = data, method = mixed())
#'
#' # Fit model using the CG method
#' mod <- fitGEV(y ~ pb(x), data = data, method = CG())
#'
#' # Fremantle data from the ismev package
#' ismev <- requireNamespace("ismev", quietly = TRUE)
#'
#' if (ismev) {
#'   data(fremantle, package = "ismev")
#'   # Transform Year
#'   fremantle <- transform(fremantle, cYear = Year - median(Year))
#'
#'   mod <- fitGEV(SeaLevel ~ pb(SOI), data = fremantle)
#'   plot(fremantle$SOI, fremantle$SeaLevel)
#'   lines(fremantle$SOI, fitted(mod))
#'
#'   mod <- fitGEV(SeaLevel ~ pb(cYear), data = fremantle)
#'   plot(fremantle$cYear, fremantle$SeaLevel)
#'   lines(fremantle$cYear, fitted(mod))
#'
#'   mod <- fitGEV(SeaLevel ~ pb(cYear) + pb(SOI), data = fremantle)
#'   plot(fremantle$cYear, fremantle$SeaLevel)
#'   lines(fremantle$cYear, fitted(mod))
#'
#'   mod <- fitGEV(SeaLevel ~ SOI, data = fremantle)
#'   plot(fremantle$SOI, fremantle$SeaLevel)
#'   lines(fremantle$SOI, fitted(mod))
#'
#'   mod <- fitGEV(SeaLevel ~ cYear, data = fremantle)
#'   plot(fremantle$cYear, fremantle$SeaLevel)
#'   lines(fremantle$cYear, fitted(mod))
#'
#'   mod <- fitGEV(SeaLevel ~ SOI + cYear, data = fremantle)
#'   plot(fremantle$SOI, fremantle$SeaLevel)
#'   lines(fremantle$SOI, fitted(mod))
#'   plot(fremantle$cYear, fremantle$SeaLevel)
#'   lines(fremantle$cYear, fitted(mod))
#' }
#' @export
fitGEV <- function(formula, data, scoring = c("fisher", "quasi"),
                   stepLength = 1, stepAttempts = 2, stepReduce = 2, eps, ...) {
  # Check that one of the correct values of scoring has been supplied
  scoring <- match.arg(scoring)
  # Choose the scoring algorithm and links
  # Fit using the supplied/default value of step length
  if (scoring == "fisher") {
    mod <- try(gamlss::gamlss(formula = formula,
                              family = GEVfisher(...),
                              mu.step = stepLength, sigma.step = stepLength,
                              nu.step = stepLength, data = data, ...),
               silent = TRUE)
  } else {
    mod <- try(gamlss::gamlss(formula = formula,
                              family = GEVquasiNewton(...),
                              mu.step = stepLength, sigma.step = stepLength,
                              nu.step = stepLength, data = data, ...),
               silent = TRUE)
  }
  # If an error is thrown then try again stepAttempts times, each  time reducing
  # the step length by a factor of a half
  isError <- inherits(mod, "try-error")
  while(isError & stepAttempts >= 1) {
    stepLength <- stepLength / stepReduce
    if (scoring == "fisher") {
      mod <- try(gamlss::gamlss(formula = formula,
                                family = GEVfisher(...),
                                mu.step = stepLength, sigma.step = stepLength,
                                nu.step = stepLength, data = data, ...),
                 silent = TRUE)
    } else {
      mod <- try(gamlss::gamlss(formula = formula,
                                family = GEVquasiNewton(...),
                                mu.step = stepLength, sigma.step = stepLength,
                                nu.step = stepLength, data = data, ...),
                 silent = TRUE)
    }
#    mod <- try(gamlss::gamlss(formula = formula, family = algor,
#                              mu.step = stepLength, sigma.step = stepLength,
#                              nu.step = stepLength, data = data, ...),
#               silent = TRUE)
    isError <- inherits(mod, "try-error")
    stepAttempts <- stepAttempts - 1L
  }
  class(mod) <- c("gamlssx", class(mod))
  return(mod)
}

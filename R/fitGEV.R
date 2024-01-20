#' Fit a Generalized extreme value (GEV) GAMLSS
#'
#' Describe  [`GEV`]
#'
#' @param formula Argument `formula` passed to [`gamlss`][`gamlss::gamlss`].
#' @param stepLength A positive numeric scalar. The initial value of the common
#'  step lengths `mu.step`, `sigma.step` and `nu.step` passed to
#'  [`gamlss.control`][`gamlss::gamlss.control`] in the first attempt to fit
#'  the model by calling [`gamlss`][`gamlss::gamlss`].
#' @param stepReduce A non-negative integer. If the first call to
#'  [`gamlss`][`gamlss::gamlss`] throws an error then we make `stepReduce`
#'  further attempts to fit the model, each time dividing by \eqn{2} the values
#'  of `mu.step`, `sigma.step` and `nu.step` supplied to
#'  [`gamlss.control`][`gamlss::gamlss.control`]. If `stepReduce < 1` then
#'  no further attempts are made.
#' @param ... Further arguments passed to [`gamlss`][`gamlss::gamlss`].
#'
#' @details Add details. Explain `stepReduce` in more detail.
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
#' @export
fitGEV <- function(formula, stepLength = 1, stepReduce = 2, ...)  {
  mod <- try(gamlss(y ~ pb(x), family = GEV, mu.step = stepLength,
                    sigma.step = stepLength, nu.step = stepLength, ...),
             silent = TRUE)
  # If an error is thrown then try again stepReduce times, each  time reducing
  # the step length by a factor of a half
  isError <- inherits(mod, "try-error")
  while(isError & stepReduce >= 1) {
    stepLength <- stepLength / 2L
    mod <- gamlss(y ~ pb(x), family = GEV, mu.step = stepLength,
                  sigma.step = stepLength, nu.step = stepLength, ...)
    isError <- inherits(mod, "try-error")
    stepReduce <- stepReduce - 1L
    stepLength <- stepLength / 2
  }
  class(mod) <- c("gamlssx", class(mod))
  return(mod)
}

#' Fit a Generalized Extreme value (GEV) GAMLSS model
#'
#' Describe
#'
#' @inheritParams gamlss::gamlss
#'
#' @param scoring A character scalar. If `scoring = "fisher"` then the weights
#'   used in the fitting algorithm are based on the expected Fisher
#'   information, that is, a Fisher's scoring algorithm is used.
#'   If `scoring = "quasi"` then these weights are based on the cross products
#'   of the first derivatives of the log-likelihood, leading to a quasi Newton
#'   scoring algorithm.
#' @param mu.link,sigma.link,xi.link Character scalars to set the respective
#'   link functions for the location (`mu`), scale (`sigma`) and shape (`xi`)
#'   parameters. The latter is passed to [`gamlss`][`gamlss::gamlss`] as
#'   `nu.link`.
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
#' @param ... Further arguments passed to [`gamlss`][`gamlss::gamlss`], in
#'   particular `method`, which sets the fitting algorithm, with options
#'   `RS()`, `CG()` or `mixed()`. The default, `method = RS()` seems to work
#'   well, as does `method = mixed()`. In contrast, `method = CG()` often
#'   requires the step length to be reduced before convergence is achieved.
#'   `fitGEV()` attempts to do this automatically. See `stepAttemmpts`.
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
#' # Transform Year
#' fremantle <- transform(fremantle, cYear = Year - median(Year))
#'
#' mod <- fitGEV(SeaLevel ~ pb(SOI), data = fremantle)
#' plot(fremantle$SOI, fremantle$SeaLevel)
#' lines(fremantle$SOI, fitted(mod))
#'
#' mod <- fitGEV(SeaLevel ~ pb(cYear), data = fremantle)
#' plot(fremantle$cYear, fremantle$SeaLevel)
#' lines(fremantle$cYear, fitted(mod))
#'
#' mod <- fitGEV(SeaLevel ~ pb(cYear) + pb(SOI), data = fremantle)
#' plot(fremantle$cYear, fremantle$SeaLevel)
#' lines(fremantle$cYear, fitted(mod))
#'
#' mod <- fitGEV(SeaLevel ~ SOI, data = fremantle)
#' plot(fremantle$SOI, fremantle$SeaLevel)
#' lines(fremantle$SOI, fitted(mod))
#'
#' mod <- fitGEV(SeaLevel ~ cYear, data = fremantle)
#' plot(fremantle$cYear, fremantle$SeaLevel)
#' lines(fremantle$cYear, fitted(mod))
#'
#' mod <- fitGEV(SeaLevel ~ SOI + cYear, data = fremantle)
#' plot(fremantle$SOI, fremantle$SeaLevel)
#' lines(fremantle$SOI, fitted(mod))
#' plot(fremantle$cYear, fremantle$SeaLevel)
#' lines(fremantle$cYear, fitted(mod))
#' @export
fitGEV <- function(formula, data, scoring = c("fisher", "quasi"),
                   mu.link = "identity", sigma.link = "log",
                   xi.link = "identity", stepLength = 1, stepAttempts = 2,
                   stepReduce = 2, ...) {
  # Check that one of the correct values of scoring has been supplied
  scoring <- match.arg(scoring)
  # Set the scoring algorithm and links
  # For all the gamlss methods to work on the returned fitted model object, we
  # need the call to gamlss::gamlss to include explicitly the names of the
  # links for mu, sigma and nu (xi here) as character scalars.
  # To achieve this, we create the internal function templateFit() and
  # modify the body of this function to include the names of the link
  # functions. This is rather clunky, but it works! Other attempts at passing
  # the link functions do not work completely. For example, vcov.gamlss(object)
  # does not work, because it performs calculations using the fitted model
  # object and needs to know the link functions.
  #
  # Fit using the supplied/default value of step length
  if (scoring == "fisher") {
    algor <- substitute(
      GEVfisher(mu.link = mu.link, sigma.link = sigma.link, nu.link = xi.link)
      )
  } else {
    algor <- substitute(
      GEVquasi(mu.link = mu.link, sigma.link = sigma.link, nu.link = xi.link)
      )
  }
  # Add the link functions to the call to gamlss() in fisherFit()
  templateFit <- function(formula, stepLength, data, ...) {
    dangerous <- NULL
    return(dangerous)
  }
  body(templateFit)[[2]] <- substitute(
    dangerous <- try(gamlss::gamlss(formula = formula, family = algor,
                                    mu.step = stepLength,
                                    sigma.step = stepLength,
                                    nu.step = stepLength, data = data, ...),
                     silent = TRUE)
    )
  cat("stepLength = ", stepLength, "\n")
  mod <- templateFit(formula = formula, stepLength = stepLength, data = data,
                     ...)
  # If an error is thrown then try again stepAttempts times, each time reducing
  # the step length by a factor of a halfs
  isError <- inherits(mod, "try-error")
  while(isError & stepAttempts >= 1) {
    stepLength <- stepLength / stepReduce
    cat("stepLength = ", stepLength, "\n")
    # We need to update the value of stepLength in templateFit()
    body(templateFit)[[2]] <- substitute(
      dangerous <- try(gamlss::gamlss(formula = formula, family = algor,
                                      mu.step = stepLength,
                                      sigma.step = stepLength,
                                      nu.step = stepLength, data = data, ...),
                       silent = TRUE)
    )
    mod <- templateFit(formula = formula, stepLength = stepLength, data = data,
                       ...)
    isError <- inherits(mod, "try-error")
    stepAttempts <- stepAttempts - 1L
  }
  if (isError) {
    stop("No convergence. An error was thrown from the last call to gamlss()")
  }
  class(mod) <- c("gamlssx", class(mod))
  return(mod)
}

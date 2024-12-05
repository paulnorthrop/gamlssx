#' Fit a Generalized Extreme value (GEV) GAMLSS Model
#'
#' Fits a Generalized Additive Model (GAM) for Location, Scale and Shape with
#' a GEV response distribution, using the function
#' [`gamlss::gamlss()`][`gamlss::gamlss`].
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
#'   parameters. The latter is passed to [`gamlss::gamlss()`][`gamlss::gamlss`]
#'   as `nu.link`.
#' @param stepLength A numeric vector of positive values. The initial
#'    values of the step lengths `mu.step`, `sigma.step` and `nu.step` passed to
#'   [`gamlss::gamlss.control()`][`gamlss::gamlss.control`] in the first attempt
#'   to fit the model by calling [`gamlss::gamlss()`][`gamlss::gamlss`]. If
#'   `stepLength` has a length that is less than 3 then `stepLength` is
#'   recycled to have length 3.
#' @param stepAttempts A non-negative integer. If the first call to
#'   [`gamlss::gamlss()`][`gamlss::gamlss`] throws an error then we make
#'   `stepAttempts` further attempts to fit the model, each time dividing by 2
#'   the values of `mu.step`, `sigma.step` and `nu.step` supplied to
#'   [`gamlss::gamlss.control()`][`gamlss::gamlss.control`]. If
#'   `stepAttempts < 1` then no further attempts are made.
#' @param stepReduce A number greater than 1. The factor by which the step
#'   lengths in `stepLength` are reduced for each extra attempt to fit the
#'   model. The default, `stepReduce = 2` means that the step lengths are
#'   halved for each extra attempt.
#' @param steps A logical scalar. Pass `steps = TRUE` to write to the
#'   console the current value of `stepLength` for each call to
#'   [`gamlss::gamlss()`][`gamlss::gamlss`].
#' @param ... Further arguments passed to
#'   [`gamlss::gamlss()`][`gamlss::gamlss`], in particular `method`, which sets
#'   the fitting algorithm, with options `RS()`, `CG()` or `mixed()`. The
#'   default, `method = RS()` seems to work well, as does `method = mixed()`.
#'   In contrast, `method = CG()` often requires the step length to be reduced
#'   before convergence is achieved. `fitGEV()` attempts to do this
#'   automatically. See `stepAttempts`. Pass `trace = FALSE`
#'   (to [`gamlss::gamlss.control()`][`gamlss::gamlss.control`]) to avoid
#'   writing to the console the global deviance after each outer iteration of
#'   the gamlss fitting algorithm.
#' @details See [`gamlss::gamlss()`][`gamlss::gamlss`] for information about
#'   the model and the fitting algorithm.
#' @return Returns a `gamlss` object. See the **Value** section of
#'   [`gamlss::gamlss()`][`gamlss::gamlss`]. The class of the returned object is
#'   `c("gamlssx", "gamlss", "gam", "glm", "lm")`.
#' @seealso [`GEV`],
#'   [`gamlss.dist::gamlss.family()`][`gamlss.dist::gamlss.family`],
#'   [`gamlss::gamlss()`][`gamlss::gamlss`]
#' @examples
#' # Load gamlss, for the function pb()
#' library(gamlss)
#'
#' ##### Simulated data
#'
#' set.seed(17012023)
#' n <- 100
#' x <- stats::runif(n)
#' mu <- 1 + 2 * x
#' sigma <- 1
#' xi <- 0.25
#' y <- nieve::rGEV(n = 1, loc = mu, scale = sigma, shape = xi)
#' data <- data.frame(y = as.numeric(y), x = x)
#' plot(x, y)
#'
#' # Fit model using the default RS method with Fisher's scoring
#' mod <- fitGEV(y ~ gamlss::pb(x), data = data)
#' # Summary of model fit
#' summary(mod)
#' # Residual diagnostic plots
#' plot(mod, xlab = "x", ylab = "y")
#' # Data plus fitted curve
#' plot(data$x, data$y, xlab = "x", ylab = "y")
#' lines(data$x, fitted(mod))
#'
#' # Fit model using the mixed method and quasi-Newton scoring
#' # Use trace = FALSE to prevent writing the global deviance to the console
#' mod <- fitGEV(y ~ pb(x), data = data, method = mixed(), scoring = "quasi",
#'               trace = FALSE)
#'
#' # Fit model using the CG method
#' # The default step length of 1 needs to be reduced to enable convergence
#' # Use steps = TRUE to write the step lengths to the console
#' mod <- fitGEV(y ~ pb(x), data = data, method = CG(), steps = TRUE)
#'
#' ##### Fremantle annual maximum sea levels
#' ##### See also the gamlssx package README file
#'
#' # Transform Year so that it is centred on 0
#' fremantle <- transform(fremantle, cYear = Year - median(Year))
#'
#' # Plot sea level against year and against SOI
#' plot(fremantle$Year, fremantle$SeaLevel, xlab = "year", ylab = "sea level (m)")
#' plot(fremantle$SOI, fremantle$SeaLevel, xlab = "SOI", ylab = "sea level (m)")
#'
#' # Fit a model with P-spline effects of cYear and SOI on location and scale
#' # The default links are identity for location and log for scale
#' mod <- fitGEV(SeaLevel ~ pb(cYear) + pb(SOI),
#'              sigma.formula = ~ pb(cYear) + pb(SOI),
#'              data = fremantle)
#'
#' # Summary of model fit
#' summary(mod)
#' # Model diagnostic plots
#' plot(mod)
#' # Worm plot
#' wp(mod)
#' # Plot of the fitted component smooth functions
#' # Note: gamlss::term.plot() does not include uncertainty about the intercept
#' # Location mu
#' term.plot(mod, rug = TRUE, pages = 1)
#' # Scale sigma
#' term.plot(mod, what = "sigma", rug = TRUE, pages = 1)
#' @export
fitGEV <- function(formula, data, scoring = c("fisher", "quasi"),
                   mu.link = "identity", sigma.link = "log",
                   xi.link = "identity", stepLength = 1, stepAttempts = 2,
                   stepReduce = 2, steps = FALSE,  ...) {
  # Check that one of the correct values of scoring has been supplied
  scoring <- match.arg(scoring)
  # Force stepLength to have length 3
  stepLength <- rep_len(stepLength, 3)
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
                                    mu.step = stepLength[1],
                                    sigma.step = stepLength[2],
                                    nu.step = stepLength[3], data = data,
                                    ...),
                     silent = TRUE)
    )
  if (steps) {
    cat("stepLength =", stepLength, "\n")
  }
  mod <- templateFit(formula = formula, stepLength = stepLength, data = data,
                     ...)
  # If an error is thrown then try again stepAttempts times, each time reducing
  # the step length by a factor of a half
  isError <- inherits(mod, "try-error")
  while(isError & stepAttempts >= 1) {
    stepLength <- stepLength / stepReduce
    if (steps) {
      cat("stepLength =", stepLength, "\n")
    }
    # We need to update the value of stepLength in templateFit()
    body(templateFit)[[2]] <- substitute(
      dangerous <- try(gamlss::gamlss(formula = formula, family = algor,
                                      mu.step = stepLength[1],
                                      sigma.step = stepLength[2],
                                      nu.step = stepLength[3], data = data,
                                      ...),
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

#' gamlssx: Generalized Additive Extreme Value Models for Location, Scale and
#' Shape
#'
#' Fits generalized additive models for the location, scale and shape
#' parameters of a generalized extreme value response distribution. The
#' methodology is based on Rigby and Stasinopoulos (2005) and implemented using
#' functions from the `gamlss` package \doi{10.32614/CRAN.package.gamlss}.
#'
#' @details The main function in `gamlssx` is [`fitGEV()`], which calls the
#'   function [`gamlss::gamlss()`][`gamlss::gamlss`].
#'   See the \href{https://paulnorthrop.github.io/gamlssx/}{gamlssx package
#'   page on Github} for more information.
#'
#' @references Rigby R.A. and Stasinopoulos D.M. (2005). Generalized additive
#'   models for location, scale and shape (with discussion), *Appl. Statist.*,
#'   **54**, part 3, pages 507-554. \doi{10.1111/j.1467-9876.2005.00510.x}
#' @seealso [`fitGEV()`], [`gamlss::gamlss()`][`gamlss::gamlss`]
#' @import gamlss
#' @importFrom stats sd
#' @docType package
"_PACKAGE"

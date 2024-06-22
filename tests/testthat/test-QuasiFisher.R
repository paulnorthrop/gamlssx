# Check that for fitGEV() scoring = "fisher" (default) and scoring = "quasi"
# produce equivalent results

tol <- 1e-3

# gamlss

# fisher

modFisher <- gamlssx::fitGEV(SeaLevel ~ 1, data = fremantle)
mFisher <- as.numeric(modFisher$mu.coefficients)
sFisher <- as.numeric(exp(modFisher$sigma.coefficients))
xFisher <- as.numeric(modFisher$nu.coefficients)
fisherEstimates <- c(mFisher, sFisher, xFisher)

# quasi

modQuasi <- gamlssx::fitGEV(SeaLevel ~ 1, data = fremantle, scoring = "quasi")
mQuasi <- as.numeric(modQuasi$mu.coefficients)
sQuasi <- as.numeric(exp(modQuasi$sigma.coefficients))
xQuasi <- as.numeric(modQuasi$nu.coefficients)
quasiEstimates <- c(mQuasi, sQuasi, xQuasi)

test_that("fisher logLik equals quasi logLik", {
  testthat::expect_equal(logLik(modFisher), logLik(modQuasi), tolerance = tol)
})

test_that("gamlss estimates equal ismev estimates", {
  testthat::expect_equal(fisherEstimates, quasiEstimates, tolerance = tol)
})

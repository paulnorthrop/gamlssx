# Check that when the default fit fails, which is usually when the CG algorithm
# has been used, the reduction in step lengths implemented in fitGEV() can help

tol <- 1e-2

# Simulate some data
set.seed(17012023)
n <- 100
x <- stats::runif(n)
mu <- 1 + 2 * x
sigma <- 1
xi <- 0.25
y <- nieve::rGEV(n = 1, loc = mu, scale = sigma, shape = xi)
plot(x, y)
data <- data.frame(y = as.numeric(y), x = x)

# Fit model using the default RS method
modRS <- fitGEV(y ~ pb(x), data = data)
muhat <- modRS$mu.coefficients
sigmahat <- exp(modRS$sigma.coefficients)
xihat <- modRS$nu.coefficients
RSestimates <- as.numeric(c(muhat, sigmahat, xihat))
# 1.2726165, 1.3830108 0.9161198 0.1926386

# Fit model using the CG method
modCG <- fitGEV(y ~ pb(x), data = data, method = CG())
muhat <- modCG$mu.coefficients
sigmahat <- exp(modCG$sigma.coefficients)
xihat <- modCG$nu.coefficients
CGestimates <- as.numeric(c(muhat, sigmahat, xihat))


test_that("RS logLik equals CG logLik", {
  testthat::expect_equal(logLik(modRS), logLik(modCG), tolerance = tol)
})

test_that("RS estimates equal CG estimates", {
  testthat::expect_equal(RSestimates, CGestimates, tolerance = tol)
})

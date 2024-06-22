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
y <- gamlssx::rGEV(n = 1, mu = mu, sigma = sigma, nu = xi)
data <- data.frame(y = as.numeric(y), x = x)

# Fit model using the default RS method
modRS <- fitGEV(y ~ gamlss::pb(x), data = data)
print(modRS)
muhat <- modRS$mu.coefficients
sigmahat <- exp(modRS$sigma.coefficients)
xihat <- modRS$nu.coefficients
RSestimates <- as.numeric(c(muhat, sigmahat, xihat))
loglikRS <- logLik(modRS)

# Fit model using the CG method
modCG <- fitGEV(y ~ gamlss::pb(x), data = data, method = CG())
muhat <- modCG$mu.coefficients
sigmahat <- exp(modCG$sigma.coefficients)
xihat <- modCG$nu.coefficients
CGestimates <- as.numeric(c(muhat, sigmahat, xihat))
loglikCG <- logLik(modCG)


test_that("RS logLik equals CG logLik", {
  testthat::expect_equal(loglikRS, loglikCG, tolerance = tol)
})

test_that("RS estimates equal CG estimates", {
  testthat::expect_equal(RSestimates, CGestimates, tolerance = tol)
})

# For this example, we need extra attempts if using the CG algorithm
# Set (extra) stepAttempts to 0 to see this

# Fit model using the CG method
test_that("CG needs extra iterations with a reduced step length", {
  testthat::expect_error(fitGEV(y ~ gamlss::pb(x), data = data, method = CG(), stepAttempts = 0))
})

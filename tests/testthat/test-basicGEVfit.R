# Check that a basic GEV fit with no covariates gives the correct results

tol <- 1e-4

# gamlss

mod <- gamlssx::fitGEV(SeaLevel ~ 1, data = fremantle, steps = TRUE)
# Estimated coefficients
muhat <- as.numeric(mod$mu.coefficients) # 1.482304
sigmahat <- as.numeric(exp(mod$sigma.coefficients)) # 0.1412218
xihat <- as.numeric(mod$nu.coefficients) # -0.2172112
# Fit again using these as starting estimates
# I pass the numeric values for mu.start etc because when the next line of
# code is run outside of test_that(), using devtools::check(), if I pass
# mu.start = muhat then the value of muhat is not passed correctly to
# templatefit(...) in fitGEV() and an error that says that muhat is not found
# is thrown. It works fine using devtools::test() and when executing code via
# the command line. I'm not entirely sure what the problem is - environments,
# presumbly. This is a quick fix.
mod <- gamlssx::fitGEV(SeaLevel ~ 1, data = fremantle, mu.start = 1.482304,
                       sigma.start = 0.1412218, nu.start = -0.2172112)
# Estimated coefficients
muhat <- mod$mu.coefficients
sigmahat <- exp(mod$sigma.coefficients)
xihat <- mod$nu.coefficients
# Estimates, SEs, maximised log-likelihood
gamlssEstimates <- as.numeric(c(muhat, sigmahat, xihat))
ses <- sqrt(diag(vcov(mod)))
ses[2] <- ses[2] * gamlssEstimates[2]
gamlssSEs <- as.numeric(ses)
gamlssLoglik <- as.numeric(logLik(mod))

# ismev

#x <- ismev::gev.fit(xdat = fremantle$SeaLevel)
#x$mle
#x$se
#-x$nllh

ismevEstimates <- c(1.4823409,  0.1412671, -0.2174320)
ismevSEs <- c(0.01672502, 0.01149461, 0.06377394)
ismevLoglik <- 43.56663

test_that("gamlss logLik equals ismev logLik", {
  testthat::expect_equal(gamlssLoglik, ismevLoglik, tolerance = tol)
})

test_that("gamlss estimates equal ismev estimates", {
  testthat::expect_equal(gamlssEstimates, ismevEstimates, tolerance = tol)
})

test_that("gamlss SEs equal ismev SEs", {
  testthat::expect_equal(gamlssSEs, ismevSEs, tolerance = tol)
})

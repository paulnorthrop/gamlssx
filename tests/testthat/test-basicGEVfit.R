# Check that a basic GEV fit with no covariates gives the correct results

tol <- 1e-4

# gamlss

mod <- fitGEV(SeaLevel ~ 1, data = fremantle)
# Estimated coefficients
muhat <- coef(mod, what = "mu")
sigmahat <- exp(coef(mod, what = "sigma"))
xihat <- coef(mod, what = "nu")
# Fit again using these as starting estimates
mod <- fitGEV(SeaLevel ~ 1, data = fremantle, mu.start = muhat,
              sigma.start = sigmahat, nu.start = xihat)
# Estimated coefficients
muhat <- coef(mod, what = "mu")
sigmahat <- exp(coef(mod, what = "sigma"))
xihat <- coef(mod, what = "nu")
# Estimates, SEs, maximised log-likelihood
gamlssEstimates <- as.numeric(c(muhat, sigmahat, xihat))
ses <- sqrt(diag(vcov(mod)))
ses[2] <- ses[2] * estimates[2]
gamlssSEs <- as.numeric(ses)
gamlssLoglik <- as.numeric(logLik(mod))

# ismev

x <- ismev::gev.fit(xdat = fremantle$SeaLevel)
x$mle
x$se
-x$nllh
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

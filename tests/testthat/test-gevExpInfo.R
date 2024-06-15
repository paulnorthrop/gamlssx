# Check that the calculation of the expected information for the GEV
# distribution agrees with mev::gev.infomat()

sigma <- 1
mevMat <- cbind(
c(1/sigma^2, -0.422784335098467/sigma^2, 0.41184033042644/sigma),
c(-0.422784335098467/sigma^2, 1.82368066085288/sigma^2, 0.332484907160274/sigma),
c(0.41184033042644/sigma, 0.332484907160274/sigma, 2.42360605517703)
)

# Based on the constants
evilsMat1 <- cbind(c(gev11e(1, 0), gev12e0Constant, gev13e0Constant),
                   c(gev12e0Constant, gev22e0Constant, gev23e0Constant),
                   c(gev13e0Constant, gev23e0Constant, gev33e0Constant))

# Based on the functions that create the constants
evilsMat1f <- cbind(c(gev11e(1, 0), gev12e0Fn(), gev13e0Fn()),
                   c(gev12e0Fn(), gev22e0Fn(), gev23e0Fn()),
                   c(gev13e0Fn(), gev23e0Fn(), gev33e0Fn()))

# Based on calling the functions with xi = 0
evilsMat2 <- cbind(c(gev11e(1, 0), gev12e(1, 0), gev13e(1, 0)),
                   c(gev12e(1, 0), gev22e(1, 0), gev23e(1, 0)),
                   c(gev13e(1, 0), gev23e(1, 0), gev33e(0)))

test_that("Expected Information constants agree with mev::gev.infomat()", {
  testthat::expect_equal(mevMat, evilsMat1)
})

test_that("Expected Information constant funs agree with mev::gev.infomat()", {
  testthat::expect_equal(mevMat, evilsMat1f)
})

test_that("Expected Information constants agree with mev::gev.infomat()", {
  testthat::expect_equal(mevMat, evilsMat2)
})

test_that("gevExpInfo(shape = 0) agrees with mev::gev.infomat()", {
  testthat::expect_equal(mevMat, gevExpInfo(scale = 1, shape = 0),
                         ignore_attr = TRUE)
})

# Check that gevExpInfo() throws an error when shape <= -0.5

test_that("gevExpInfo() errors for shape = -1/2", {
  testthat::expect_error(gevExpInfo(scale = 1, shape = -1/2))
})
test_that("gevExpInfo() errors for shape = -1", {
  testthat::expect_error(gevExpInfo(scale = 1, shape = -1))
})

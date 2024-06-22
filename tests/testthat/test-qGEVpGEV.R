# Check that pGEV() and qGEV() are consistent

## Example shape parameters

# Positive
xi1 <- 0.1
low <- -1 / xi1

# Zero
xi2 <- 0

# Negative
xi3 <- -1e-7
up <- -1 / xi3

## Example input vector

pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)

test_that("pgev and qgev are consistent", {
  expect_equal(pGEV(qGEV(pvec, nu = xi1), nu = xi1), pvec)
  expect_equal(pGEV(qGEV(pvec, nu = xi2), nu = xi2), pvec)
  expect_equal(pGEV(qGEV(pvec, nu = xi3), nu = xi3), pvec)
})

context("Test all the components of package")

## estimate_kevin is correct

test_that("estimate_kevin works", {
  set.seed(10)
  y <- rnorm(10)
  res <- estimate_kevin(y)

  expect_true(length(res) == 1)
  expect_true(res >= 1)
  expect_true(res < length(y))
  expect_true(res %% 1 == 0)
})

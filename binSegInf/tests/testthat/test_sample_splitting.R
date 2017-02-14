context("Test sample splitting")

## sample_splitting is correct

test_that("sample_splitting works with binSeg_fixedSteps", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 50)) + 0.1*rnorm(100)
  
  res <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  
  expect_true(class(res) == "ss")
  expect_true(res$method == "binSeg_fixedSteps")
  expect_true(res$jumps == 49)
})

test_that("sample_splitting works with fLasso_fixedSteps", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 50)) + 0.1*rnorm(100)
  
  res <- sample_splitting(y, fLasso_fixedSteps, numSteps = 1)
  
  expect_true(class(res) == "ss")
  expect_true(res$method == "fLasso_fixedSteps")
  expect_true(res$jumps == 49)
})
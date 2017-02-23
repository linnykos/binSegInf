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

###################################

## jumps.ss is correct

test_that("jumps.ss works", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 50)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  res <- jumps(obj)
  
  expect_true(res == 49)
})

##################################

## contrast_vector_ss is correct

test_that("contrast_vector_ss works", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 50)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  
  res <- contrast_vector_ss(obj, 1)
  expect_true(all(which(res < 0) == seq(2, 48, by = 2)))
  expect_true(all(which(res > 0) == seq(52, 100, by = 2)))
})

#####################################

## pvalue_ss is correct

test_that("pvalue_ss works", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 50)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  v <- contrast_vector_ss(obj, 1)
  
  res <- pvalue_ss(y, v)
  expect_true(length(res) == 1)
  expect_true(res < 0.05)
})
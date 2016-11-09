context("Test binary segmentation with fixed step size")

## .cusum_contrast is correct

test_that(".cusum_contrast returns a vector", {
  res <- .cusum_contrast(1, 5, 10)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == 10)
})

test_that(".cusum_contrast returns a vector with the right values", {
  res <- .cusum_contrast(1, 5, 10)
  
  expect_true(sum(res) == 0)
  expect_true(length(unique(res)) == 2)
  expect_true(length(unique(res[1:5])) == 1)
  expect_true(length(unique(res[6:10])) == 1)
})

######################

## .cusum is correct

test_that(".cusum returns 0 for the flat vector", {
  res <- .cusum(rep(0, 10), 1, 5, 10)
  expect_true(res == 0)
})

#######################

## .find_breakpoint is correct

test_that(".find_breakpoint is correct", {
  y <- c(rep(0,5), rep(1,5))
  res <- .find_breakpoint(y, 1, 10)
  
  expect_true(length(res) == 2)
  expect_true(res$breakpoint == 5)
})

test_that(".find_breakpoint reports 0 if start = end", {
  y <- c(rep(0,5), rep(1,5))
  res <- .find_breakpoint(y, 4, 4)
  
  expect_true(length(res) == 2)
  expect_true(res$breakpoint == 4)
  expect_true(res$cusum == 0)
})

#############################

## binSeg_fixedSteps is correct

test_that("binSeg_fixedSteps works on one jump", {
  y <- c(rep(0, 10), rep(1, 10))
  res <- binSeg_fixedSteps(y, 1)
  
  expect_true(length(res) == 2)
  expect_true(class(res) == "bsFs")
  expect_true(class(res$tree)[1] == "Node")
  expect_true(res$tree$breakpoint == 10)
})

test_that("binSeg_fixedSteps works with two jumps", {
  y <- c(rep(0, 10), rep(5, 5), rep(6, 5))
  res <- binSeg_fixedSteps(y, 2)
  
  expect_true(res$tree$breakpoint == 10)
  expect_true(res$tree$children[[2]]$breakpoint == 15)
})

test_that("the isValid function works", {
  y <- c(rep(0, 10), rep(1, 10))
  res <- binSeg_fixedSteps(y, 1)
  
  expect_true(isValid(res))
})

#################################

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

test_that(".cusum returns negatives", {
  res <- .cusum(c(rep(0,5), rep(-5, 5)), 1, 5, 10)
  
  expect_true(res < 0)
})

#######################

## .find_breakpoint is correct

test_that(".find_breakpoint is correct", {
  set.seed(10)
  y <- c(rep(0,5), rep(1,5)) + 0.01*rnorm(10)
  res <- .find_breakpoint(y, 1, 10)
  
  expect_true(length(res) == 2)
  expect_true(res$breakpoint == 5)
})

test_that(".find_breakpoint reports 0 if start = end", {
  set.seed(10)
  y <- c(rep(0,5), rep(1,5)) + 0.01*rnorm(10)
  res <- .find_breakpoint(y, 4, 4)
  
  expect_true(length(res) == 2)
  expect_true(res$breakpoint == 4)
  expect_true(res$cusum == 0)
})

test_that(".find_breakpoint splits at the best location", {
  set.seed(10)
  y <- c(rep(0, 5), rep(10,4), -9) + 0.01*rnorm(10)
  
  mat <- cbind(1, 1:9, 10)
  cusum.vec <- apply(mat, 1, function(x){
    .cusum(y, x[1], x[2], x[3])
  })
  
  res <- .find_breakpoint(y, 1, 10)
  
  idx <- which.max(abs(cusum.vec))
  
  expect_true(cusum.vec[idx] == res$cusum)
})

test_that(".find_breakpoint reports negatives", {
  set.seed(10)
  y <- c(rep(0, 5), rep(-10, 5)) + 0.01*rnorm(10)
  res <- .find_breakpoint(y, 1, 10)
   
  expect_true(res$cusum < 0)
})

#############################

## binSeg_fixedSteps is correct

test_that("binSeg_fixedSteps works on one jump", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  res <- binSeg_fixedSteps(y, 1)
  
  expect_true(length(res) == 2)
  expect_true(class(res) == "bsFs")
  expect_true(class(res$tree)[1] == "Node")
  expect_true(res$tree$breakpoint == 10)
})

test_that("binSeg_fixedSteps works with two jumps", {
  set.seed(10)
  y <- c(rep(0, 10), rep(5, 5), rep(6, 5)) + 0.01*rnorm(20)
  res <- binSeg_fixedSteps(y, 2)
  
  expect_true(res$tree$breakpoint == 10)
  expect_true(res$tree$children[[2]]$breakpoint == 15)
})

test_that("the is_valid function works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  res <- binSeg_fixedSteps(y, 1)
  
  expect_true(is_valid(res))
})

test_that("binSeg errors on duplicated values", {
  y <- rep(0, 10)
  expect_error(binSeg_fixedSteps(y, 1))
})

test_that("binSeg splits at the best location", {
  set.seed(10)
  y <- c(rep(0, 5), rep(10,4), -9) + 0.01*rnorm(10)
  obj <- binSeg_fixedSteps(y, 1)
  
  breakpoint <- obj$tree$breakpoint
  
  mat <- cbind(1, 1:9, 10)
  cusum.vec <- apply(mat, 1, function(x){
    .cusum(y, x[1], x[2], x[3])
  })
  
  expect_true(all(abs(cusum.vec[breakpoint]) >= abs(cusum.vec)))
})

test_that("binSeg works with many jumps", {
  set.seed(10)
  y <- rnorm(100)
  res <- binSeg_fixedSteps(y, 10)
  
  expect_true(length(res$tree$leaves) == 11)
})

#################################

## .cusum_contrast_full is correct

test_that(".cusum_contrast_full works", {
  res <- .cusum_contrast_full(3, 7, 10, 15)
  
  expect_true(length(res) == 15)
  expect_true(all(res[-c(3:10)] == 0))
  expect_true(all(res[3:10] != 0))
})

############################

## jumps.bsFs is correct

test_that("jumps.bsFs works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(5, 5), rep(6, 5)) + 0.01*rnorm(20)
  obj <- binSeg_fixedSteps(y, 2)
  
  res <- jumps(obj)
  
  expect_true(all(res == c(10, 15)))
})

#################################

## jumps_cusum.bsFs is correct

test_that("jumps_cusum.bsFs works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(5, 5), rep(6, 5)) + 0.01*rnorm(20)
  obj <- binSeg_fixedSteps(y, 2)
  
  res <- jump_cusum(obj)
  
  expect_true(length(res) == 2)
  expect_true(all(res > 0))
})

test_that("jumps_cusum.bsFs reports negative values",{
  set.seed(10)
  y <- c(rep(0, 10), rep(-5, 5), rep(-10, 5)) + 0.01*rnorm(20)
  obj <- binSeg_fixedSteps(y, 2)
  
  res <- jump_cusum(obj)
  
  expect_true(all(res < 0))
})
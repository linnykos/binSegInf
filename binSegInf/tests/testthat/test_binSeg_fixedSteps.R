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

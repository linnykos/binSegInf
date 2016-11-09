context("Test binary segmentation with fixed step size")

## .cusum_contrast is correct

test_that(".cusum_contrast returns a vector", {
  res <- .cusum_contrast(1, 5, 10)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == 10)
})

test_that(".cusum_contrast returns a vector with the right values", {
  res <- .cusum_contrast(1, 5, 10)
  
  expect_true(all(res == c(rep(-1/5,5), rep(1/5,5))))
})

######################

## .cusum is correct

test_that(".cusum returns 0 for the flat vector", {
  res <- .cusum(rep(0, 10), 1, 5, 10)
  expect_true(res == 0)
})
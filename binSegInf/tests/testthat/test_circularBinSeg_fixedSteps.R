context("Test circular binary segmentation fixed steps")

## .cusum_cbs is correct

test_that(".cusum_cbs works", {
  vec <- cumsum(1:10)
  x <- c(3,5)
  res <- .cusum_cbs(x, vec)
  
  answer <- sqrt(1/3 + 1/7) * (mean(3:5) - mean(c(1:2,6:10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when the left side is length 1", {
  vec <- cumsum(1:10)
  x <- c(2,5)
  res <- .cusum_cbs(x, vec)
  
  answer <- sqrt(1/4 + 1/6) * (mean(2:5) - mean(c(1,6:10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when the right side is length 1", {
  vec <- cumsum(1:10)
  x <- c(6,9)
  res <- .cusum_cbs(x, vec)
  
  answer <- sqrt(1/4 + 1/6) * (mean(6:9) - mean(c(1:5,10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when the right side is length 1", {
  vec <- cumsum(1:10)
  x <- c(6,9)
  res <- .cusum_cbs(x, vec)
  
  answer <- sqrt(1/4 + 1/6) * (mean(6:9) - mean(c(1:5,10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when there is no left side", {
  vec <- cumsum(1:10)
  x <- c(1,5)
  res <- .cusum_cbs(x, vec)
  
  answer <- sqrt(1/5 + 1/5) * (mean(1:5) - mean(c(6:10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when there is no right side", {
  vec <- cumsum(1:10)
  x <- c(7,10)
  res <- .cusum_cbs(x, vec)
  
  answer <- sqrt(1/4 + 1/6) * (mean(7:10) - mean(c(1:6)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs fails when the bounds are the entire interval", {
  vec <- cumsum(1:10)
  x <- c(1,10)
  expect_error(.cusum_cbs(x, vec))
})

######################


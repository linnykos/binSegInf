context("Test circular binary segmentation fixed steps")

## .cusum_cbs is correct

test_that(".cusum_cbs works", {
  vec <- cumsum(1:10)
  x <- c(3,5)
  res <- .cusum_cbs(x, vec)
  
  answer <- sqrt(1/(1/3 + 1/7)) * (mean(3:5) - mean(c(1:2,6:10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when the left side is length 1", {
  vec <- cumsum(1:10)
  x <- c(2,5)
  res <- .cusum_cbs(x, vec)
  
  answer <- sqrt(1/(1/4 + 1/6)) * (mean(2:5) - mean(c(1,6:10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when the right side is length 1", {
  vec <- cumsum(1:10)
  x <- c(6,9)
  res <- .cusum_cbs(x, vec)
  
  answer <- sqrt(1/(1/4 + 1/6)) * (mean(6:9) - mean(c(1:5,10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when the right side is length 1", {
  vec <- cumsum(1:10)
  x <- c(6,9)
  res <- .cusum_cbs(x, vec)
  
  answer <- sqrt(1/(1/4 + 1/6)) * (mean(6:9) - mean(c(1:5,10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when there is no left side", {
  vec <- cumsum(1:10)
  x <- c(1,5)
  res <- .cusum_cbs(x, vec)
  
  answer <- sqrt(1/(1/5 + 1/5)) * (mean(1:5) - mean(c(6:10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when there is no right side", {
  vec <- cumsum(1:10)
  x <- c(7,10)
  res <- .cusum_cbs(x, vec)
  
  answer <- sqrt(1/(1/4 + 1/6)) * (mean(7:10) - mean(c(1:6)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs fails when the bounds are the entire interval", {
  vec <- cumsum(1:10)
  x <- c(1,10)
  expect_error(.cusum_cbs(x, vec))
})

test_that(".cusum_cbs computes correctly when there are negative values", {
  set.seed(10)
  y <- c(rnorm(10), rnorm(10, mean = 10), rnorm(10))
  vec <- cumsum(y)
  
  x <- c(11, 11)
  res <- .cusum_cbs(x, vec)
  answer <- sqrt(1/(1/1 + 1/29))*mean(y[11] - y[c(1:10, 12:30)])
  expect_true(abs(res - answer) <= 1e-6)
  
  x <- c(11, 20)
  res <- .cusum_cbs(x, vec)
  answer <- sqrt(1/(1/10 + 1/20))*mean(y[11:20] - y[c(1:10, 21:30)])
  expect_true(abs(res - answer) <= 1e-6)
})

######################

## .enumerate_breakpoints_cbs is correct

test_that(".enumerate_breakpoints_cbs works", {
  n <- 5
  res <- .enumerate_breakpoints_cbs(n)
  expect_true(is.matrix(res))
  expect_true(ncol(res) == 2)
  expect_true(all(res[,1] <= res[,2]))
  
  map <- res[,1]*6+res[,2]
  expect_true(length(map) == length(unique(map))) #all unique rows
})

test_that(".enumerate_breakpoints_cbs outputs the right rows", {
  n <- 4
  res <- .enumerate_breakpoints_cbs(n)
  answer <- matrix(c(1,1, 1,2, 1,3, 2,2, 2,3, 2,4, 3,3, 3,4, 4,4), ncol = 2, byrow = T)
  
  expect_true(nrow(res) == nrow(answer))
  
  map_res <- res[,1]*5+res[,2]
  map_answer <- answer[,1]*5+answer[,2]
  expect_true(all(sort(map_res) == sort(map_answer)))
})

####################

## .find_breakpoint_cbs is correct

test_that(".find_breakpoint_cbs works", {
  vec <- cumsum(1:10)
  res <- .find_breakpoint_cbs(vec)
  
  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(length(res$breakpoint) == 2)
  expect_true(res$breakpoint[1] >= 1)
  expect_true(res$breakpoint[2] <= 10)
  expect_true(res$breakpoint[1] <= res$breakpoint[2])
  expect_true(is.numeric(res$cusum))
})

test_that(".find_breakpoint_cbs finds the right breakpoint", {
  set.seed(10)
  y <- c(rnorm(10), rnorm(10, mean = 10), rnorm(10))
  vec <- cumsum(y)
  res <- .find_breakpoint_cbs(vec)
  
  expect_true(all(res$breakpoint == c(11,20)))
})

####################

## circularBinSeg_fixedSteps is correct

test_that("circularBinSeg_fixedSteps works", {
  set.seed(10)
  y <- c(rnorm(10), rnorm(10, mean = 10), rnorm(10))
  res <- circularBinSeg_fixedSteps(y, 1)
  
  expect_true(class(res) == "cbsFs")
})

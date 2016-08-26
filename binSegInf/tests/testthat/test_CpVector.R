context("Test changepoint")

## test .computeJumpIdx

test_that("it splits the jump indices evenly", {
  n <- 100
  jump.loc <- seq(0,1,length.out = 5)[-c(1,5)]
  res <- .computeJumpIdx(n, jump.loc)
  d <- diff(diff(res))
  
  expect_true(all(d <= 1))
  expect_true(all(d >= -1))
})

##################

## test .formMeanVec

test_that("it forms the mean vector properly", {
  jump.loc <- c(25, 50, 75)
  jump.height <- 1:4
  res <- .formMeanVec(100, jump.loc, jump.height)
  
  expect_true(length(res) == 100)
  
  d <- diff(table(res))
  expect_true(all(d <= 1))
  expect_true(all(d >= -1))
})

######################

## test CpVector

test_that("it forms a proper CpVector class", {
  res <- CpVector(n = 100, 1:4, c(.25,.5,.75))
  expect_true(length(res$data) == 100)
  expect_true(class(res) == "CpVector")
  
  expect_true(all(res$jump.height == 1:4))
  expect_true(all(res$jump.idx == c(25,50,75)))
})
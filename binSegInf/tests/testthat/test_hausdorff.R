context("Test Hausdorff")

## test hausdorff

test_that("hausdorff works properly", {
  set1 <- c(1,5,10)
  set2 <- set1 + 1
  
  expect_true(hausdorff(set1, set2) == 1)
})

#########################

## test enumerateJumps

test_that("jumps works properly", {
  vec <- rep(1:4, each = 10)
  res <- enumerateJumps(vec)
  
  expect_true(all(res == c(10,20,30)))
})

test_that("it works when vec is length 1", {
  vec <- 5
  expect_true(length(enumerateJumps(vec)) == 0)
})

test_that("it works when vec is all singleton", {
  vec <- 1:10
  res <- enumerateJumps(vec)
  
  expect_true(all(res == 1:9))
})
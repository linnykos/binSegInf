context("Test Hausdorff")

## test hausdorff

test_that("it works properly", {
  set1 <- c(1,5,10)
  set2 <- set1 + 1
  
  expect_true(hausdorff(set1, set2) == 1)
})
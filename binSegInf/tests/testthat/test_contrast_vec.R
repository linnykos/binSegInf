context("Test contrast vector")

## contrast_vector.bsFs is correct

test_that("contrast_vector.bsFs works", {
  y <- c(rep(0, 10), rep(5, 5), rep(6, 5))
  obj <- binSeg_fixedSteps(y, 2)
  
  res <- contrast_vector(obj, 1)
  expect_true(length(res) == 20)
  expect_true(sum(res) == 0)
  expect_true(all(res[1:10] == -1/10))
  expect_true(all(res[11:15] == 1/5))
})
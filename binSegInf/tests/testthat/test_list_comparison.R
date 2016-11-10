context("Test list_comparison")

test_that(".list_comparison works", {
  y <- c(rep(0, 20), rep(5, 10), rep(10, 5), rep(11, 5))
  obj <- binSeg_fixedSteps(y, 3)
  
  res <- .list_comparison(obj)
  
  expect_true(length(res) == 3)
  expect_true(is.list(res))
})
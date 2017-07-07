context("Test circular binSeg fixed Step polyhedra")

## polyhedra.cbsFs is correct

test_that("polyhedra.cbsFs works", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedSteps(y, 2)
  res <- polyhedra(obj)
  
  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(sort(names(res)) == c("gamma", "u")))
  expect_true(ncol(res$gamma) == length(y))
  expect_true(length(res$u) == nrow(res$gamma))
})

test_that("polyhedra.cbsFs gives the right number of rows", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedSteps(y, 2)
  res <- polyhedra(obj)
  
  n <- length(y)
  jump_vec <- jumps(obj, sorted = F)
  start <- jump_vec[1]; end <- jump_vec[2]
  
  ans <- 2*((nrow(.enumerate_breakpoints_cbs(n))-1) + 
    (nrow(.enumerate_breakpoints_cbs(start)) + nrow(.enumerate_breakpoints_cbs(n-end)) +
    nrow(.enumerate_breakpoints_cbs(end - start)) - 1))
  
  expect_true(nrow(res$gamma) == ans)
})
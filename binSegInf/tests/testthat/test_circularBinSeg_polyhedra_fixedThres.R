context("Test circular binSeg fixed Threshold polyhedra")

## .get_signs_cbsFt is correct

test_that(".get_signs_cbsFt works", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedThres(y, thres = 10)
  res <- .get_signs_cbsFt(obj)
  
  expect_true(length(res) == 7)
  expect_true(all(res[3:7] == 0))
  expect_true(all(res[1:2] != 0))
})

############

## polyhedra.cbsFt is correct

test_that("polyhedra.cbsFt works", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedThres(y, thres = 10)
  res <- polyhedra(obj)
  
  expect_true(class(res) == "polyhedra")
  expect_true(all(names(res) == c("gamma", "u")))
  expect_true(ncol(res$gamma) == length(y))
  expect_true(nrow(res$gamma) == length(res$u))
})

test_that("polyhedra.cbsFt satisfies the polyhedra requirement", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedThres(y, thres = 10)
  res <- polyhedra(obj)

  expect_true(all(res$gamma %*% y >= res$u))
})

test_that("polyhedra.cbsFt satisfies the polyhedra requirement", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedThres(y, thres = 10)
  res <- polyhedra(obj)

  expect_true(all(res$gamma %*% y >= res$u))
})
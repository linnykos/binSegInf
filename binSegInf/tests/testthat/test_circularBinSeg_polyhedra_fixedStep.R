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
  expect_true(all(res$u == 0))
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

test_that("polyhedra.cbsFs satisfies polyhedra requirement", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedSteps(y, 2)
  res <- polyhedra(obj)
  
  expect_true(all(res$gamma %*% y >= res$u))
})


test_that("having the same model if and only if the inequalities are satisfied", {
  set.seed(5)
  y <- c(rep(0,10), rep(-2,10), rep(-1,5)) + rnorm(25)
  obj <- circularBinSeg_fixedSteps(y,2)
  
  model <- jumps(obj, sorted = F)
  poly <- polyhedra(obj)
  
  expect_true(all(poly$gamma %*% y >= poly$u))
  
  trials <- 100
  for(i in 1:trials){
    set.seed(i*10)
    y_tmp <- c(rep(0,10), rep(-2,10), rep(-1,5)) + rnorm(25)
    obj_tmp <- circularBinSeg_fixedSteps(y_tmp,2)
    
    model_tmp <- jumps(obj_tmp, sorted = F)
    
    bool1 <- all(length(model) == length(model_tmp) && all(model == model_tmp))
    bool2 <- all(poly$gamma %*% y_tmp >= poly$u)
    
    expect_true(bool1 == bool2)
  }
})

test_that("having the same model if and only if the inequalities are satisfied", {
  set.seed(5)
  y <- rnorm(25)
  obj <- circularBinSeg_fixedSteps(y,2)
  
  model <- obj$model[,1:2]
  poly <- polyhedra(obj)
  
  expect_true(all(poly$gamma %*% y >= poly$u))
  
  trials <- 100
  for(i in 1:trials){
    set.seed(i*10)
    y_tmp <- c(rnorm(5), rnorm(5, mean = -10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
    obj_tmp <- circularBinSeg_fixedSteps(y_tmp,2)
    
    model_tmp <- jumps(obj_tmp, sorted = F)
    
    bool1 <- all(length(model) == length(model_tmp) && all(model == model_tmp))
    bool2 <- all(poly$gamma %*% y_tmp >= poly$u)
    
    expect_true(bool1 == bool2)
  }
})

######################

## .cusum_cbs_contrast_full is correct

test_that(".cusum_cbs_contrast_full works", {
  res <- .cusum_cbs_contrast_full(1, c(5,10), 15, 15)
  
  expect_true(length(res) == 15)
})

test_that(".cusum_cbs_contrast_full gives the correct answer", {
  set.seed(10)
  y <- rnorm(15)
  con <- .cusum_cbs_contrast_full(1, c(5,10), 15, 15)
  cusum <- .cusum_cbs(c(5,10), cumsum(y))
  
  expect_true(abs(con%*%y - cusum) <= 1e-6)
})

test_that(".cusum_cbs_contrast_full gives the correct answer with shift", {
  set.seed(10)
  y <- rnorm(15)
  con <- .cusum_cbs_contrast_full(5, c(7,9), 10, 15)
  cusum <- .cusum_cbs(c(3,5), cumsum(y[5:10]))
  
  expect_true(abs(con%*%y - cusum) <= 1e-6)
})

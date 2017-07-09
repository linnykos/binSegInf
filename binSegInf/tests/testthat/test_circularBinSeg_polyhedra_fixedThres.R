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

test_that("polyhedra.cbsFt satisfies the polyhedra requirement for 1 jump", {
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


test_that("polyhedra.cbsFt having the same model if and only if the inequalities are satisfied", {
  set.seed(5)
  y <- c(rep(0,10), rep(-2,10), rep(-1,5)) + rnorm(25)
  obj <- circularBinSeg_fixedThres(y,2)
  
  model <- jumps(obj, sorted = F)
  poly <- polyhedra(obj)
  
  expect_true(all(poly$gamma %*% y >= poly$u))
  
  trials <- 100
  for(i in 1:trials){
    set.seed(i*10)
    y_tmp <- c(rep(0,10), rep(-2,10), rep(-1,5)) + rnorm(25)
    obj_tmp <- circularBinSeg_fixedThres(y_tmp,2)
    
    model_tmp <- jumps(obj_tmp, sorted = F)
    
    bool1 <- all(length(model) == length(model_tmp) && all(model == model_tmp))
    bool2 <- all(poly$gamma %*% y_tmp >= poly$u)
    
    expect_true(bool1 == bool2)
  }
})

test_that("polyhedra.cbsFt having the same model if and only if the inequalities are satisfied, wrong model", {
  set.seed(5)
  y <- rnorm(25)
  obj <- circularBinSeg_fixedThres(y,2)
  
  model <- obj$model[,1:2]
  poly <- polyhedra(obj)
  
  expect_true(all(poly$gamma %*% y >= poly$u))
  
  trials <- 100
  for(i in 1:trials){
    set.seed(i*10)
    y_tmp <- c(rnorm(5), rnorm(5, mean = -10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
    obj_tmp <- circularBinSeg_fixedThres(y_tmp,2)
    
    model_tmp <- jumps(obj_tmp, sorted = F)
    
    bool1 <- all(length(model) == length(model_tmp) && all(model == model_tmp))
    bool2 <- all(poly$gamma %*% y_tmp >= poly$u)
    
    expect_true(bool1 == bool2)
  }
})
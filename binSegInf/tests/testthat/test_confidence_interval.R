context("Test confidence interval")

## confidence_interval is correct

test_that("confidence_interval has the right coverage zero", {
  set.seed(10)
  y <- rnorm(20)
  obj <- binSeg_fixedSteps(y, 1)
  
  poly <- polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res <- confidence_interval(y, poly, contrast, gridsize = 100)
  expect_true(0 >= res[1]-.5)
  expect_true(0 <= res[2]+.5)
})

test_that("confidence_interval has right coverage non-zero", {
  set.seed(10)
  y <-  c(rep(0, 10), rep(3, 10)) + rnorm(20)
  obj <- binSeg_fixedSteps(y, 1)
  
  poly <- polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res <- confidence_interval(y, poly, contrast, gridsize = 100)
  expect_true(3 >= res[1]-.5)
  expect_true(3 <= res[2]+.5)
})

test_that("confidence interval one and two-sided are related", {
  set.seed(10)
  y <-  c(rep(0, 10), rep(3, 10)) + rnorm(20)
  obj <- binSeg_fixedSteps(y, 1)
  
  poly <- polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res.twosided <- confidence_interval(y, poly, contrast, gridsize = 50)
  res.pos.onesided <- confidence_interval(y, poly, contrast, gridsize = 50, 
    alternative = "one.sided")
  res.neg.onesided <- confidence_interval(y, poly, -contrast, gridsize = 50, 
    alternative = "one.sided")
  
  expect_true(res.twosided[1] <= res.pos.onesided[1])
  expect_true(res.pos.onesided[1] >= res.neg.onesided[1])
  expect_true(res.neg.onesided[1] <= 0.5)
  expect_true(res.pos.onesided[2] == Inf)
  expect_true(res.neg.onesided[2] == Inf)
})

test_that("confidence interval is not a point", {
  set.seed(14)
  
  dat <- CpVector(100, 0, NA)
  y <- dat$data
  
  obj <- binSeg_fixedSteps(y, 1)

  poly <- polyhedra(obj, y)
  contrast <- contrast_vector(obj, 1)

  res <- confidence_interval(y, poly, contrast, gridsize = 100)
  
  expect_true(abs(res[1]-res[2]) > 1e-4)
})

test_that("confidence int. should give a left point less than right", {
  set.seed(1)
  dat <- CpVector(100, 0, NA)
  y <- dat$data
  
  obj <- binSeg_fixedSteps(y, 1)

  poly <- polyhedra(obj, y)
  contrast <- contrast_vector(obj, 1)

  res <- confidence_interval(y, poly, contrast, gridsize = 50)
  
  expect_true(res[1] <= res[2])
})

test_that("confindence int throws warning appropriate", {
  set.seed(500)
  dat <- CpVector(100, 0, NA)
  y <- dat$data
  
  obj <- binSeg_fixedSteps(y, 1)

  poly <- polyhedra(obj, y)
  contrast <- contrast_vector(obj, 1)

  expect_warning(confidence_interval(y, poly, contrast, gridsize = 50))
})

#####################################

## .select_index is correct

test_that(".select_index selects the correct left index", {
  vec <- dnorm(seq(-3, 3, length.out = 50))
  vec[40] <- 0.05
  
  res <- .select_index(vec, 0.05)
  expect_true(res != 40)
  expect_true(all(vec[1:res] <= 0.05))
})

test_that(".select_index selects the correct right index", {
  vec <- pnorm(seq(-3, 3, length.out = 50))
  vec[30] <- 0.95
  
  res <- .select_index(vec, 0.95, F)
  expect_true(res != 30)
  expect_true(all(vec[res:50] >= 0.95))
})

test_that(".select_index throws warning if vec and alpha don't overlap", {
  vec <- 1:10
  alpha <- 20
  expect_warning(.select_index(vec, alpha))
})
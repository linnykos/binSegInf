context("Test confidence interval")

## confidence_interval is correct

test_that("confidence_interval has the right coverage zero", {
  set.seed(10)
  y <- rnorm(20)
  obj <- binSeg_fixedSteps(y, 1)
  
  poly <- form_polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res <- confidence_interval(y, poly, contrast, gridsize = 100)
  expect_true(0 >= res[1]-.5)
  expect_true(0 <= res[2]+.5)
})

test_that("confidence_interval has right coverage non-zero", {
  set.seed(10)
  y <-  c(rep(0, 10), rep(3, 10)) + rnorm(20)
  obj <- binSeg_fixedSteps(y, 1)
  
  poly <- form_polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res <- confidence_interval(y, poly, contrast, gridsize = 100)
  expect_true(3 >= res[1]-.5)
  expect_true(3 <= res[2]+.5)
})

test_that("confidence interval one and two-sided are related", {
  set.seed(10)
  y <-  c(rep(0, 10), rep(3, 10)) + rnorm(20)
  obj <- binSeg_fixedSteps(y, 1)
  
  poly <- form_polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res.twosided <- confidence_interval(y, poly, contrast, gridsize = 50)
  res.pos.onesided <- confidence_interval(y, poly, contrast, gridsize = 50, 
    alternative = "one.sided")
  res.neg.onesided <- confidence_interval(y, poly, -contrast, gridsize = 50, 
    alternative = "one.sided")
  
  expect_true(res.twosided[1] >= res.pos.onesided[1])
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

  poly <- form_polyhedra(obj, y)
  contrast <- contrast_vector(obj, 1)

  res <- confidence_interval(y, poly, contrast, gridsize = 100)
  
  expect_true(abs(res[1]-res[2]) > 1e-4)
})
context("Test confidence interval")

## confidence_interval is correct

test_that("confidence_interval has the right coverage zero", {
  set.seed(10)
  y <- rnorm(20)
  obj <- binSeg_fixedSteps(y, 1)
  
  poly <- form_polyhedra(obj, y)
  contrast <- contrast_vector(obj, 1)
  
  res <- confidence_interval(y, poly, contrast, gridsize = 100)
  expect_true(0 >= res[1]-.5)
  expect_true(0 <= res[2]+.5)
})

test_that("confidence_interval has right coverage non-zero", {
  set.seed(10)
  y <-  c(rep(0, 10), rep(3, 10)) + rnorm(20)
  obj <- binSeg_fixedSteps(y, 1)
  
  poly <- form_polyhedra(obj, y)
  contrast <- contrast_vector(obj, 1)
  
  res <- confidence_interval(y, poly, contrast, gridsize = 100)
  expect_true(3 >= res[1]-.5)
  expect_true(3 <= res[2]+.5)
})
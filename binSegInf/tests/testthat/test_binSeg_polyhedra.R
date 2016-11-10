context("Test binSeg_polyhedra")

## form_polyhedra.bsFs is correct

test_that("form_polyhedra.bsFs works", {
  y <- c(rep(0, 10), rep(10, 5), rep(5, 5))
  obj <- binSeg_fixedSteps(y, 2)
})

###############################

## .vector_matrix_signedDiff is correct

test_that(".vector_matrix_signedDiff works", {
  res <- .vector_matrix_signedDiff(c(1,2,3), matrix(1:6, ncol = 3, byrow = T),
    1, c(1,-1))
  
  expect_true(all(res == matrix(c(0,0,0,5,7,9), ncol = 3, nrow = 2, byrow = T)))
})

test_that(".vector_matrix_signedDiff gets signs right", {
  res <- .vector_matrix_signedDiff(rep(0,3), matrix(1:15, ncol = 3), 
    0, c(1,1,1,-1,-1))
  
  expect_true(all(apply(res, 1, function(x){length(unique(sign(x)))}) == 1))
})

test_that(".vector_matrix_signedDiff returns a matrix of the right dim", {
  res <- .vector_matrix_signedDiff(1:10, matrix(1:60, ncol = 10), 
    1, rep(c(-1,1), each = 3))
  
  expect_true(all(dim(res) == c(6, 10)))
})

###################################

## .gammaRows_from_comparisons is correct

test_that(".gammaRows_from_comparisons works", {
  set.seed(10)
  y <- sample(c(1:10))
  vec <- matrix(c(1,5,10), ncol = 3)
  mat <- cbind(1, c(1:9)[-5], 10)
  
  res <- .gammaRows_from_comparisons(vec, mat, 10, y)
  
  expect_true(all(dim(res) == c(17, 10)))
})

test_that(".gammaRow_from_comparisons is fulfilled by y", {
  y <- c(rep(0, 5), rep(10,4), -9)
  obj <- binSeg_fixedSteps(y, 1)
  
  expect_true(obj$tree$breakpoint == 9)
  
  vec <- matrix(c(1,9,10), ncol = 3)
  mat <- cbind(1, 1:8, 10)
  
  res <- .gammaRows_from_comparisons(vec, mat, 10, y)
  
  expect_true(all(res %*% y >= 0))
})
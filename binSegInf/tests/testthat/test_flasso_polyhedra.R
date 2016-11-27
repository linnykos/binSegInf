context("Test fused lasso polyhedra")

## .compute_fused_numerator_polyhedra is correct

test_that(".compute_fused_numerator_polyhedra returns a matrix of correct size", {
  D <- .form_Dmatrix(10)
  res <- .compute_fused_numerator_polyhedra(D, c(1:4,6:7,9))
  
  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == c(9-2, 10)))
})

###################################

## .form_contrast_flasso is correct

test_that(".form_contrast_flasso returns the right matrices", {
  mat <- matrix(1:30, 6, 5)
  vec <- 1:6
  sign.win <- 1
  active.idx <- 2
  
  res <- .form_contrast_flasso(mat, vec, sign.win, active.idx)
  
  expect_true(length(res) == 2)
  expect_true(sum(abs(res - mat[2,]/(1+2))) < 1e-7)
  expect_true(all(names(res) == c("win", "lose")))
  expect_true(is.numeric(res$win))
  expect_true(length(res$win) == 5)
  expect_true(is.matrix(res$lose))
  expect_true(all(dim(res$lose) == c(2*5, 5)))
  expect_true(!any(is.infinite(res$lose)))
})
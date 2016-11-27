context("Test fused lasso polyhedra")

## .compute_fused_numerator_polyhedra is correct

test_that(".compute_fused_numerator_polyhedra returns a matrix of correct size", {
  D <- .form_Dmatrix(10)
  res <- .compute_fused_numerator_polyhedra(D, c(1:4,6:7,9))
  
  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == c(9-2, 10)))
})
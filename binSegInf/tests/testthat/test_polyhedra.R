context("Test polyhedra")

## polyhedra is correct

test_that("polyhedra forms a correct class", {
  res <- polyhedra(matrix(1:25,5,5), rep(1,5))
  
  expect_true(class(res) == "polyhedra")
})

########################

## isValid.polyhedra is correct

test_that("isValid.polyhedra works", {
  res <- polyhedra(matrix(1:25,5,5), rep(1,5))
  
  expect_true(isValid(res))
  
  res$u <- rep(1,4)
  
  expect_error(isValid(res))
})
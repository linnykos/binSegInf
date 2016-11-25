context("Test polyhedra")

## polyhedra is correct

test_that("polyhedra forms a correct class", {
  res <- polyhedra(matrix(1:25,5,5), rep(1,5))
  
  expect_true(class(res) == "polyhedra")
})

########################

## is_valid.polyhedra is correct

test_that("is_valid.polyhedra works", {
  res <- polyhedra(matrix(1:25,5,5), rep(1,5))
  
  expect_true(is_valid(res))
  
  res$u <- rep(1,4)
  
  expect_error(is_valid(res))
})
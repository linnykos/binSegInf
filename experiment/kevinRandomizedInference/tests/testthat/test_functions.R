context("Test all the components of package")

## estimate_kevin is correct

test_that("estimate_kevin works", {
  set.seed(10)
  y <- rnorm(10)
  res <- estimate_kevin(y)

  expect_true(length(res) == 1)
  expect_true(res >= 1)
  expect_true(res < length(y))
  expect_true(res %% 1 == 0)
})

##################

## .polyhedron_vector_generator is correct

test_that(".polyhedron_vector_generator works", {
  res <- .polyhedron_vector_generator(4, 10)

  expect_true(all(res == c(rep(1/4, 4), rep(-1/6, 6))))
})

#################

## polyhedron_kevin is correct

test_that("polyhedron_kevin works", {
  y <- rep(0, 10)
  i <- 7
  res <- polyhedron_kevin(y,i)

  n <- length(y)
  expect_true(all(dim(res) == c(n-2,n)))
  expect_true(is.matrix(res))
})

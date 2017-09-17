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
  expect_true(is.numeric(res))
  expect_true(!is.na(res))
})

test_that("estimate_kevin works for extreme signals", {
  set.seed(10)
  y <- c(rnorm(4, mean = 50), rnorm(6))
  res <- estimate_kevin(y)

  expect_true(res == 4)
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

test_that("polyhedron_kevin satisfies polyhedron inequality", {
  trials <- 100
  res <- sapply(1:trials, function(x){
    set.seed(x)
    y <- rnorm(10)
    i <- estimate_kevin(y)
    gam <- polyhedron_kevin(y,i)

    all(gam %*% y >= 0)
  })

  expect_true(all(res))
})

test_that("polyhedron_kevin generates Gamma specific to selection", {
  trials <- 100
  set.seed(1)
  y <- rnorm(10)
  i <- estimate_kevin(y)
  gam <- polyhedron_kevin(y,i)

  res <- sapply(1:trials, function(x){
    set.seed(10*x)
    y2 <- rnorm(10)
    i2 <- estimate_kevin(y2)

    bool1 <- i == i2
    bool2 <- all(gam%*%y2 >= 0)

    bool1 == bool2
  })

  expect_true(all(res))
})

###################



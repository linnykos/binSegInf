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

## poly.pval_kevin is correct

test_that("poly.pval_kevin works", {
  set.seed(10)
  y <- rnorm(10)
  i <- estimate_kevin(y)
  mat <- polyhedron_kevin(y, i)
  contrast <- .polyhedron_vector_generator(i, length(y))

  res <- poly.pval_kevin(mat, y, 1, contrast)

  expect_true(length(res) == 3)
  expect_true(all(is.numeric(unlist(res))))
  expect_true(res$pvalue <= 1)
  expect_true(res$pvalue >= 0)
  expect_true(all(unlist(res) >= 0))
})

test_that("poly.pval_kevin forms uniform pvalues", {
  trials <- 500
  vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    y <- rnorm(10)
    i <- estimate_kevin(y)
    mat <- polyhedron_kevin(y,i)
    contrast <- .polyhedron_vector_generator(i, length(y))

    poly.pval_kevin(mat, y, 1, contrast)$pvalue
  })

  expect_true(sum(abs(quantile(vec) - c(0, 0.25, 0.5, 0.75, 1))) <= 0.1)
})

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
  trials <- 2000
  vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    y <- rnorm(10)
    i <- estimate_kevin(y)
    mat <- polyhedron_kevin(y,i)
    contrast <- .polyhedron_vector_generator(i, length(y))

    val <- poly.pval_kevin(mat, y, 1, contrast)$pvalue
    val
  })

  # plot(sort(vec), seq(0,1,length.out = length(vec)))
  # lines(c(0,1), c(0,1), col = "red", lwd = 2)

  expect_true(sum(abs(quantile(vec) - c(0, 0.25, 0.5, 0.75, 1)))/trials <= 1e-5)
})

test_that("poly.pval_kevin does not crash when contrast is chosen independent of Gamma matrix", {
  set.seed(10)
  y <- rnorm(10)
  contrast <- c(rep(1/5,5), rep(-1/5,5))

  n <- length(y); unique_seed <- sum(abs(y))
  set.seed(2*10*unique_seed)
  y_boot <- y + stats::rnorm(n, sd = 1)
  i <- estimate_kevin(y_boot)
  mat <- polyhedron_kevin(y_boot, i)
  res <- poly.pval_kevin(mat, y_boot, 1, contrast)

  expect_true(length(res) == 3)
})

####################

## .compute_truncGaus_terms is correct

test_that(".compute_truncGaus_terms works", {
  set.seed(10)
  y <- rnorm(10)
  i <- estimate_kevin(y)
  mat <- polyhedron_kevin(y,i)
  contrast <- .polyhedron_vector_generator(i, length(y))

  res <- .compute_truncGaus_terms(y, mat, contrast, 1)

  expect_true(is.list(res))
  expect_true(length(res) == 4)
  expect_true(res$a <= res$b)
})

test_that(".compute_truncGaus_terms works when contrast is chosen independent of Gamma matrix", {
  set.seed(10)
  y <- rnorm(10)
  contrast <- c(rep(1/5,5), rep(-1/5,5))

  n <- length(y); unique_seed <- sum(abs(y))
  set.seed(2*10*unique_seed)
  y_boot <- y + stats::rnorm(n, sd = 1)
  i <- estimate_kevin(y_boot)
  mat <- polyhedron_kevin(y_boot, i)
  res <- .compute_truncGaus_terms(y_boot, mat, contrast, 1)

  expect_true(res$a <= res$b)
})

#################

## .truncated_gauss_cdf is correct

test_that(".truncated_gauss_cdf works", {
  set.seed(10)
  y <- rnorm(10)
  i <- estimate_kevin(y)
  mat <- polyhedron_kevin(y, i)
  contrast <- .polyhedron_vector_generator(i, length(y))
  terms <- .compute_truncGaus_terms(y, mat, contrast, sigma = 1)

  res <- .truncated_gauss_cdf(terms$term, terms$sigma, terms$a, terms$b)

  expect_true(length(res) == 3)
  expect_true(all(is.numeric(unlist(res))))
  expect_true(res$pvalue <= 1)
  expect_true(res$pvalue >= 0)
  expect_true(all(unlist(res) >= 0))
})

test_that(".truncated_gauss_cdf will not return 0 in this test case", {
  set.seed(1470)
  y <- rnorm(10)
  i <- estimate_kevin(y)
  mat <- polyhedron_kevin(y,i)
  contrast <- .polyhedron_vector_generator(i, length(y))
  terms <- .compute_truncGaus_terms(y, mat, contrast, sigma = 1)

  res <- .truncated_gauss_cdf(terms$term, terms$sigma, terms$a, terms$b)

  bool1 <- terms$term >= terms$a
  bool2 <- res$pvalue != 0
  expect_true(bool1 == bool2)
})

#################

## sampler_kevin is correct

test_that("sampler_kevin works", {
  set.seed(10)
  y <- rnorm(10)
  contrast <- c(rep(1/5,5), rep(-1/5,5))
  res <- sampler_kevin(y, 1, 1, 100, contrast)

  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
})

test_that("poly.pval_kevin forms uniform pvalues", {
  trials <- 1000
  contrast <- c(rep(1/5,5), rep(-1/5,5))
  doMC::registerDoMC(cores = 3)

  func <- function(i){
    print(i)
    set.seed(10*i)
    y <- rnorm(10)
    sampler_kevin(y, 1, 1, 500, contrast)
  }

  vec <- foreach::"%dopar%"(foreach::foreach(i = 1:trials),
                           func(i))

  # plot(sort(vec), seq(0,1,length.out = length(vec)))
  # lines(c(0,1), c(0,1), col = "red", lwd = 2)


  expect_true(sum(abs(quantile(vec) - c(0, 0.25, 0.5, 0.75, 1)))/trials <= 1e-5)
})

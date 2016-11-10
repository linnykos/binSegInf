context("Test list_comparison")

## .comparison_mat is correct

test_that(".comparison_mat enforces winning and losing requirements", {
  expect_error(.comparison_mat(1:2, matrix(1:15, ncol = 3, nrow = 5)))
  expect_error(.comparison_mat(1:3, matrix(1:15, ncol = 5, nrow = 3)))
})

test_that(".comparison_mat works", {
  res <- .comparison_mat(1:3, matrix(1:15, ncol = 3))
  
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("winning", "losing")))
  expect_true(all(res$winning == c(1:3)))
  expect_true(all(res$losing == matrix(1:15, ncol = 3)))
})

test_that(".comparison_mat enforces ordering on losing", {
  expect_error(.comparison_mat(1:3, matrix(15:1, ncol = 3, nrow = 5)))
  
  mat <-  matrix(1:15, ncol = 3, nrow = 5)
  mat2 <- mat
  mat2[2,1] <- 8
  expect_error(.comparison_mat(1:3, mat2))
  
  mat3 <- mat
  mat3[3,3] <- 4
  expect_error(.comparison_mat(1:3, mat3))
})

####################################

## .extract_startEnd is correct

test_that(".extract_startEnd works", {
  res <- .extract_startEnd("1-10")
  expect_true(all(res == c(1,10)))
})
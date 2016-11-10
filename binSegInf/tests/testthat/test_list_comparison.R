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

#####################################

## .get_leaves_matrix_excluding is correct

test_that(".get_leaves_matrix_excluding works", {
  y <- c(rep(0, 20), rep(5, 10), rep(10, 5), rep(11, 5))
  obj <- binSeg_fixedSteps(y, 3)
  res <- .get_leaves_matrix_excluding(obj$tree, "1-40")
  
  expect_true(all(res == matrix(c(21, 40, 31, 40), 2, 2)))
})

test_that(".get_leaves_matrix_excluding errors if node is not a split", {
  y <- c(rep(0, 20), rep(5, 10), rep(10, 5), rep(11, 5))
  obj <- binSeg_fixedSteps(y, 3)
  
  expect_error(.get_leaves_matrix_excluding(obj$tree, "1-23"))
})

#########################################

## .threeColumnMatrix_from_nodeVec is correct

test_that(".threeColumnMatrix_from_nodeVec works", {
  res <- .threeColumnMatrix_from_nodeVec(c(1,5))
  
  expect_true(all(res == cbind(1, 1:4, 5)))
})

test_that(".threeColumnMatrix_from_nodeVec can exclude a value", {
  res <- .threeColumnMatrix_from_nodeVec(c(1,5), 4)
  
  expect_true(all(dim(res) == c(3,3)))
  expect_true(all(res == cbind(1, 1:3, 5)))
})

############################################

## .threeColumnMatrix_from_nodeMatrix is correct

test_that(".threeColumnMatrix_from_nodeMatrix works", {
  mat <- matrix(1:10, ncol = 5, nrow = 2, byrow = T)
  res <- .threeColumnMatrix_from_nodeMatrix(mat)
  
  expect_true(length(res) == 5)
  expect_true(is.list(res))
  expect_true(all(sapply(res, ncol) == 3))
  expect_true(all(sapply(res, nrow) == 5))
  
  for(i in 1:5){
    expect_true(all(res[[i]] == cbind(i, i:(i+4), i+5)))
  }
})
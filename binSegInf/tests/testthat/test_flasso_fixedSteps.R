context("Test fused lasso fixed steps")

## flasso_fixedSteps is correct

test_that("flasso_fixedSteps returns a valid matrix for 2 jumps", {
  set.seed(10)
  y <- rnorm(20)
  res <- flasso_fixedSteps(y, 2)
  
  expect_true(class(res) == "flFs")
  expect_true(all(dim(res$model) == c(2,3)))
  expect_true(all(colnames(res$model) == c("Index", "Sign", "Lambda")))
  expect_true(!any(is.na(res$model)))
  expect_true(res$numSteps == 2)
})

test_that("flasso finds the right jump", {
  set.seed(10)
  y <- c(rep(0, 10), rep(10,10)) + 0.01*rnorm(20)
  res <- flasso_fixedSteps(y, 1)
  
  expect_true(res$model[,"Index"] == 10)
  expect_true(res$model[,"Sign"] == 1)
})

test_that("flasso works with many jumps", {
  set.seed(10)
  y <- rnorm(100)
  res <- flasso_fixedSteps(y, 10)
  
  expect_true(all(dim(res$model) == c(10,3)))
  expect_true(all(res$model$Lambda == sort(res$model$Lambda, decreasing = T)))
})

test_that("flasso gest the right 2 jumps", {
  set.seed(10)
  y <- c(rep(0, 10), rep(10, 10), rep(0, 10)) + 0.01*rnorm(30)
  res <- flasso_fixedSteps(y, 2)
  
  idx <- order(res$model$Index)
  expect_true(all(res$model$Index[idx] == c(10,20)))
  expect_true(all(res$model$Sign[idx] == c(1,-1)))
  expect_true(all(res$model$Lambda == sort(res$model$Lambda, decreasing = T)))
})

test_that("flasso gets the right 4 jumps", {
  set.seed(10)
  n <- 10; h <- 50
  y <- c(rep(0,n), rep(h,n), rep(0,n), rep(h,n), rep(0,n)) + 0.01*rnorm(5*n)
  res <- flasso_fixedSteps(y, 4)
  
  expect_true(all(sort(res$model$Index) == c(10, 20, 30,40)))
  expect_true(all(res$model$Lambda == sort(res$model$Lambda, decreasing = T)))
})

##############################

## jumps.flFs is correct

test_that("jumps.flFs returns the right vector", {
  set.seed(10)
  n <- 10; h <- 50
  y <- c(rep(0,n), rep(h,n), rep(0,n), rep(h,n), rep(0,n)) + 0.01*rnorm(5*n)
  res <- flasso_fixedSteps(y, 4)
  
  expect_true(all(sort(jumps(res)) == c(10, 20, 30,40)))
})

###############################

## jump_lambda.flFs is correct

test_that("jump_lambda.flFs returns the right vector", {
  set.seed(10)
  n <- 10; h <- 50
  y <- c(rep(0,n), rep(h,n), rep(0,n), rep(h,n), rep(0,n)) + 0.01*rnorm(5*n)
  res <- flasso_fixedSteps(y, 4)
  
  expect_true(all(sort(jump_lambda(res), decreasing = T) == jump_lambda(res)))
})


#######################################

## .form_Dmatrix is correct

test_that(".form_Dmatrix forms a correct matrix", {
  res <- .form_Dmatrix(10)
  
  expect_true(all(dim(res) == c(9,10)))
  
  res2 <- matrix(0, 9, 10)
  for(i in 1:9){
    res2[i,c(i,i+1)] <- c(-1,1)
  }
  
  expect_true(all(res == res2))
})

#############################

## .select_nonactive is correct

test_that(".select_nonactive selects indices", {
  vec <- rep(NA, 10)
  vec[1:3] <- c(5,2,3)
  
  res <- .select_nonactive(20, vec)
  expect_true(length(res) == 19-3)
  expect_true(all(res == sort(res)))
  expect_true(!any(c(5,2,3) %in% res))
  expect_true(!any(duplicated(res)))
})

test_that(".select_nonactive returns full vector", {
  vec <- rep(NA, 10)
  res <- .select_nonactive(20, vec)
  expect_true(all(res == 1:19))
})

###############################

## .svd_solve is correct

test_that(".svd_solve solves correctly", {
  A <- matrix(1:25, 5, 5)
  A <- (t(A) + A)/2
  x <- c(1:5)
  b <- A%*%x
  
  res <- .svd_solve(A, b)
  
  expect_true(length(res) == 5)
  expect_true(sum(abs(x-res)) < 1e-7)
  expect_true(!is.matrix(res))
})

####################################

## .compute_fused_numerator is correct

test_that(".compute_fused_numerator returns a vector", {
  set.seed(10)
  D <- .form_Dmatrix(10)
  idx <- c(1:9)[-c(5,8)]
  y <- rnorm(10)
  
  res <- .compute_fused_numerator(D, idx, y)
  expect_true(is.numeric(res))
  expect_true(length(res) == 9 - 2)
})

####################################

## .compute_fused_denominator is correct

test_that(".compute_fused_denominator returns a vector", {
  D <- .form_Dmatrix(10)
  idx <- c(1:9)[-c(5,8)]
  model.mat <- as.data.frame(matrix(NA, 2, 2))
  model.mat[,1] <- c(5,8)
  model.mat[,2] <- c(1,-1)
  colnames(model.mat) <- c("Index", "Sign")
  
  res <- .compute_fused_denominator(D, idx, model.mat)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == 9 - 2)
})

test_that(".compute_fused_denominator can return all 0's", {
  D <- .form_Dmatrix(10)
  idx <- c(1:9)
  model.mat <- matrix(NA, 0, 2)
  colnames(model.mat) <- c("Index", "Sign")
  
  res <- .compute_fused_denominator(D, idx, model.mat)
  expect_true(length(res) == 9)
  expect_true(all(res == 0))
   
})


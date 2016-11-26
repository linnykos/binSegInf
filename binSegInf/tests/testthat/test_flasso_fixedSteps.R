context("Test fused lasso fixed steps")

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
  expect_true(length(res) == 17)
  expect_true(all(res == sort(res)))
  expect_true(!any(c(5,2,3) %in% res))
  expect_true(!any(duplicated(res)))
})

test_that(".select_nonactive returns full vector", {
  vec <- rep(NA, 10)
  res <- .select_nonactive(20, vec)
  expect_true(all(res == 1:20))
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
})
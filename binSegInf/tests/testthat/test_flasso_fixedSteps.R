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
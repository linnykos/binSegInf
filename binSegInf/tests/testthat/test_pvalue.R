context("Test pvalues from truncated Gaussians")

## pvalue is correct

test_that("p value is high power for correct changepoint", {
  set.seed(10)
  y <- c(rep(0, 10), rep(50, 10)) + rnorm(20)
  obj <- binSeg_fixedSteps(y, 1)
  
  poly <- form_polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res <- pvalue(y, poly, contrast)
  expect_true(res < .1)
})

test_that("p value are roughly uniform", {
  len <- 50
  pvalue_null.vec <- numeric(len)
  pvalue_alt.vec <- numeric(len)
  
  for(i in 1:len){
    set.seed(i*10)
    y <- rnorm(20)
    obj <- binSeg_fixedSteps(y, 1)
    
    poly <- form_polyhedra(obj)
    contrast <- contrast_vector(obj, 1)
    
    pvalue_null.vec[i] <- pvalue(y, poly, contrast)
  }
  
  for(i in 1:len){
    set.seed(i*10)
    y <-  c(rep(0, 10), rep(3, 10)) + rnorm(20)
    obj <- binSeg_fixedSteps(y, 1)
    
    poly <- form_polyhedra(obj)
    contrast <- contrast_vector(obj, 1)
    
    pvalue_alt.vec[i] <- pvalue(y, poly, contrast)
  }
  
  quant <- c(0, 0.25, 0.5, 0.75, 1)
  expect_true(sum(abs(quantile(pvalue_null.vec, probs = quant, na.rm = T) - quant))
    < sum(abs(quantile(pvalue_alt.vec, probs = quant, na.rm = T) - quant)))
})

############################

## .truncated_gauss_cdf is correct

test_that(".truncated_gauss_cdf does not give Nan", {
  res <- .truncated_gauss_cdf(10, 0, 1, 9.8, Inf)
  expect_true(res == 0)
})

###################################

## .compute_truncGaus_terms is correct
# 
# test_that(".compute_truncGaus_terms preserves vlo correctly", {
#   set.seed(10)
#   y <- c(rep(0,5), rep(-2,2), rep(-1,3)) + rnorm(10)
#   obj <- binSeg_fixedSteps(y,2)
#   
#   poly <- form_polyhedra(obj)
#   contrast <- contrast_vector(obj, 1)
# })
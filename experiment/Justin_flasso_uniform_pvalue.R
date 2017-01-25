rm(list=ls())
library(testthat)
source("tests/testthat/test_flasso_vsGenLassoInf.R")

pval_vec <- numeric(1000)
pval_vec2 <- numeric(1000)

for(i in 1:1000){
  set.seed(10*i)
  n <- 100
  y <- rnorm(n, 0, 1)
  
  #use my code
  obj <- fLasso_fixedSteps(y, 1)
  poly <- polyhedra(obj, y)
  contrast <- contrast_vector(obj, 1)
  pval_vec2[i] <- pvalue(y, poly, contrast)
  
  #use justin's code
  D <- makeDmat(n, ord = 0)
  mypath <- dualpathSvd2(y, D, 1, approx = T)
  G <- getGammat.naive(obj = mypath, y = y, condition.step = 1)
  #d <- getdvec(mypath, y, 1, 1, type = "segment")
  #pval_vec[i] <- pval.fl1d(y, G$Gammat, d, 1)
  pval_vec[i] <- pval.fl1d(y, G$Gammat, contrast, 1)
  
  if(i %% 100 == 0) cat('*')
}

par(mfrow = c(1,3))
qqplot(pval_vec, seq(0, 1, length.out = length(pval_vec)), 
       pch = 16, xlab = "Actual quantile", ylab = "Theoretical quantile")
lines(x = c(0,1), y = c(0, 1), col = "red", lwd = 2)

qqplot(pval_vec2, seq(0, 1, length.out = length(pval_vec2)), 
       pch = 16, xlab = "Actual quantile", ylab = "Theoretical quantile")
lines(x = c(0,1), y = c(0, 1), col = "red", lwd = 2)

plot(pval_vec, pval_vec2, xlab = "Justin", ylab = "Mine")

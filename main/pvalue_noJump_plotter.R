rm(list=ls())
load("../results/pvalue_noJump_bsFs_2016-11-17.RData")

res <- bsFs_0JumpPValue
n <- 100
alpha <- 0.05
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

samp.selector <- function(lis, type = NA, func = function(x){x}){
  if(is.na(type)) return(1:ncol(lis[[1]]))

  bool.vec <- sapply(lis[[1]][2,], func)
  which(bool.vec)
}

idx <- samp.selector(res[1])
hist(res[[1]][1,idx], main = "P-values", col = "gray")

qqplot(res[[1]][1,], seq(0, 1, length.out = ncol(res[[1]])), pch = 16)

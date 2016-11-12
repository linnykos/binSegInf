rm(list=ls())
load("../main/res/CI_noJump_bsFs_2016-11-11.RData")

res <- bsFs_0JumpCI
n <- 100
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

samp.selector <- function(lis, type = NA, func = function(x){x}){
  if(is.na(type)) return(1:ncol(lis[[1]]))

  bool.vec <- sapply(lis[[1]][2,], func)
  which(bool.vec)
}

idx <- samp.selector(res[1])
hist(res[[1]][1,idx], main = "P-values", col = "gray")

idx <- samp.selector(res[1], type = T, func = function(x){if(abs(x-50)<=25) TRUE else FALSE})
hist(res[[1]][1,idx], main = "P-values", col = "gray")

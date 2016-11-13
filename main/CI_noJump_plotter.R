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

coverage.vec <- apply(res[[1]][,idx], 2, function(x){
  if(x[1] >= x[2] & x[1] <= x[3]) TRUE else FALSE
})
sum(coverage.vec)/ncol(res[[1]][,idx])

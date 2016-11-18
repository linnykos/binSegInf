rm(list=ls())
load("../results/CI_noJump_bsFs_2016-11-18.RData")

res <- bsFs_0JumpCI
n <- 100
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

samp.selector <- function(lis, type = NA, func = function(x){TRUE}){
  if(is.na(type)) return(1:ncol(lis[[1]]))

  bool.vec <- sapply(lis[[1]][2,], func)
  which(bool.vec)
}

idx <- samp.selector(res[1])

coverage.vec <- apply(res[[1]][,idx], 2, function(x){
  if(x[5] != 0) return(NA)
  if(0 >= x[2] & 0 <= x[3]) TRUE else FALSE
})

coverage.vec <- coverage.vec[!is.na(coverage.vec)]
sum(coverage.vec)/length(coverage.vec)

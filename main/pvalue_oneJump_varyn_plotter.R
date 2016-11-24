rm(list=ls())
load("../results/pvalue_oneJump_varyn_bsFs_2016-11-24.RData")

res <- bsFs_1JumpPValue_nvary
vec <- numeric(length(res))
alpha <- 0.05

samp.selector <- function(lis, type = NA, func = function(x,y){x == y}){
  if(is.na(type)) return(1:ncol(lis[[1]]))
  
  jump.percent  <-  as.numeric(strsplit(names(lis), 
  split = "-")[[1]][3])
  true.loc  <- round(n*jump.percent)
  
  bool.vec <- sapply(lis[[1]][2,], func, y = true.loc)
  which(bool.vec)
}

for (i in 1:length(res)){
  idx <- samp.selector(res[i])
  vec[i] <- length(which(res[[i]][1,idx] <= alpha))/ncol(res[[i]][,idx])
}

plot(x = paramMat[,4], y = vec, pch = 16, cex = 2)
lines(x = paramMat[,4], y = vec)

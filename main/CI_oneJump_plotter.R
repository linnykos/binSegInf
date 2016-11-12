rm(list=ls())
load("../main/res/CI_oneJump_bsFs_2016-11-11.RData")

res <- bsFs_1JumpCI
n <- 100
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

mat  <- matrix (0, ncol = length(jump.loc), nrow = length (jump.height))
colnames(mat) <- as.character(jump.loc)
rownames(mat) <- as.character(jump.height)

samp.selector <- function(lis, type = NA, func = function(x){x}){
  if(is.na(type)) return(1:ncol(lis[[1]]))

  bool.vec <- sapply(lis[[1]][2,], func)
  which(bool.vec)
}

for (i in 1:length(res)){
  idx <- samp.selector(res[i])
  mat[i] <- length(intersect(which(res[[i]][2,idx] >= res[[i]][3,idx]),
    which(res[[i]][2,idx] <= res[[i]][4,idx])))/ncol(res[[i]][,idx])
}


image(.rotate(1-mat), zlim = c(0,1))

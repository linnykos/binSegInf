rm(list=ls())
load("../results/CI_oneJump_bsFs_2016-11-18.RData")

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

.detangle_matrix <- function(mat){
  x = as.numeric(colnames(mat))
  y = as.numeric(rownames(mat))
  z = t(mat)
  
  list(x = x, y = y, z = z)
}

for (i in 1:length(res)){
  idx <- samp.selector(res[i])
  mat[i] <- length(intersect(which(res[[i]][2,idx] >= res[[i]][3,idx]),
    which(res[[i]][2,idx] <= res[[i]][4,idx])))/ncol(res[[i]][,idx])
}

#check coverage
image(.detangle_matrix(mat), zlim = c(0,1))

#plot the unconditional length in absolute scale
for (i in 1:length(res)){
  idx <- samp.selector(res[i])
  na.idx <- which(res[[i]][6,] == 0)
  idx <- intersect(idx, na.idx)
  len.vec <- apply(res[[i]][,idx], 2, function(x){x[4] - x[3]})
  mat[i] <- mean(len.vec)
}

image(.detangle_matrix(mat))

#plot the unconditional length in log scale
mat2 <- mat
rownames(mat2) <- as.character(log(as.numeric(rownames(mat2))))
image(.detangle_matrix(mat2))





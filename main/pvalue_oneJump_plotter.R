rm(list=ls())
load("../results/pvalue_oneJump_bsFs_2016-11-18.RData")

res <- bsFs_1JumpPValue
n <- 100
alpha <- 0.05

mat  <- matrix (0, ncol = length(jump.loc), nrow = length (jump.height))
colnames(mat) <- as.character(jump.loc)
rownames(mat) <- as.character(jump.height)

.detangle_matrix <- function(mat){
  x = as.numeric(colnames(mat))
  y = as.numeric(rownames(mat))
  z = t(mat)
  
  list(x = x, y = y, z = z)
}

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
  mat[i] <- length(which(res[[i]][1,idx] <= alpha))/ncol(res[[i]][,idx])
}

# unconditional plot with y-axis in absolute
image(.detangle_matrix(mat), zlim = c(0,1))
contour(.detangle_matrix(mat), add = T, levels = 0.95, lwd = 3)

# unconditional plot with y-axis in log-scale
mat2 <- mat
rownames(mat2) <- as.character(log(as.numeric(rownames(mat2))))
image(.detangle_matrix(mat2), zlim = c(0,1))
contour(.detangle_matrix(mat2), add = T, levels = 0.95, lwd = 3)


## conditional on correct jump
for (i in 1:length(res)){
  idx <- samp.selector(res[i], type = T)
  mat[i] <- length(which(res[[i]][1,idx] <= alpha))/ncol(res[[i]][,idx])
}


image(.detangle_matrix(mat), zlim = c(0,1), xlim = c(0,100))
contour(.detangle_matrix(mat), add = T, levels = 0.95, lwd = 3)

rm(list=ls())
load("../results/pvalue_oneJump_varyn_2016-11-30.RData")

alpha <- 0.05

.detangle_matrix <- function(mat){
  x = as.numeric(colnames(mat))
  y = as.numeric(rownames(mat))
  z = t(mat)
  
  list(x = x, y = y, z = z)
}

samp.selector <- function(lis, true.loc, type = NA, func = function(x,y){x == y}){
  if(is.na(type)) return(1:ncol(lis[[1]]))
  
  bool.vec <- sapply(lis[[1]][2,], func, y = true.loc)
  which(bool.vec)
}

form.matrix <- function(res, alpha, resp.func = resp.power, ...){
  mat  <- matrix (0, nrow = length(jump.loc), ncol = length(n.vec))
  rownames(mat) <- as.character(jump.loc)
  colnames(mat) <- as.character(n.vec)

  for (i in 1:length(res)){
    jump.percent  <-  as.numeric(strsplit(names(res)[i],split = "-")[[1]][3])
    n <- as.numeric(strsplit(names(res)[i], split = "-")[[1]][4])
    true.loc  <- round(n*jump.percent)
    
    idx <- samp.selector(res[i], true.loc, ...)
    mat[i] <- resp.func(res[[i]], idx, true.loc, n)
  }
  
  t(mat)
}

#some default resp.func
resp.power <- function(res, idx, true.loc, n){
  length(which(res[1,idx] <= alpha))/length(idx)
}
resp.exact <- function(res, idx, true.loc, n){
  length(idx)/ncol(res)
}
resp.variance <- function(res, idx, true.loc, n){
  mean(abs(res[2,idx]-true.loc))/n
}

################################
# unconditional plot with y-axis in absolute

bsMat <- form.matrix(bsFs_1JumpPValue_nvary, alpha)
flMat <- form.matrix(flFs_1JumpPValue_nvary, alpha)

par(mfrow = c(1,2))
image(.detangle_matrix(bsMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Size of n", main = "Binary Segmentation")
contour(.detangle_matrix(bsMat), add = T, levels = 0.95, lwd = 3)

image(.detangle_matrix(flMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Size of n", main = "Fused Lasso")
contour(.detangle_matrix(flMat), add = T, levels = 0.95, lwd = 3)

########################################


#more lenient metric to getting the right jump
bsMat <- form.matrix(bsFs_1JumpPValue_nvary, alpha, resp.func = resp.variance)
flMat <- form.matrix(flFs_1JumpPValue_nvary, alpha, resp.func = resp.variance)
zlim = c(min(bsMat, flMat), max(bsMat, flMat))

rownames(bsMat) <- as.character(log(as.numeric(rownames(bsMat))))
rownames(flMat) <- as.character(log(as.numeric(rownames(flMat))))

par(mfrow = c(1,2))
image(.detangle_matrix(bsMat), zlim = zlim, xlab = "Jump Index", 
  ylab = "log Size of n", main = "Binary Segmentation")
contour(.detangle_matrix(bsMat), add = T, levels = 0.15, lwd = 3)

image(.detangle_matrix(flMat), zlim = zlim, xlab = "Jump Index", 
  ylab = "log Size of n", main = "Fused Lasso")
contour(.detangle_matrix(flMat), add = T, levels = 0.15, lwd = 3)

rm(list=ls())
load("../results/CI_oneJump_varyn_2016-12-01.RData")

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

form.matrix <- function(res, alpha, resp.func = resp.coverage, ...){
  mat  <- matrix (0, nrow = length(jump.loc), ncol = length(n.vec))
  rownames(mat) <- as.character(jump.loc)
  colnames(mat) <- as.character(n.vec)

  for (i in 1:length(res)){
    jump.percent  <-  as.numeric(strsplit(names(res)[i],split = "-")[[1]][3])
    n <- as.numeric(strsplit(names(res)[i], split = "-")[[1]][4])
    true.loc  <- round(n*jump.percent)
    
    idx1 <- samp.selector(res[i], true.loc, ...)
    idx2 <- which(res[[i]][7,] == 0)
    idx <- intersect(idx1, idx2)
    
    mat[i] <- resp.func(res[[i]], idx, true.loc, n)
  }
  
  t(mat)
}

#some default resp.func
resp.coverage <- function(res, idx, true.loc, n){
  length(intersect(which(res[2,idx] >= res[3,idx]), 
    which(res[2,idx] <= res[4,idx])))/length(idx)
}
resp.length <- function(res, idx, true.loc, n){
  len.vec <- apply(res[,idx], 2, function(x){x[4] - x[3]}); mean(len.vec)
}
resp.power <- function(res, idx, true.loc, n){
  length(unique(c(which(0 <= res[3,idx]), which(0 >= res[4,idx]))))/length(idx)
}

################################
# unconditional plot with y-axis in absolute to check coverage

bsMat <- form.matrix(bsFs_1JumpCI)
flMat <- form.matrix(flFs_1JumpCI)

par(mfrow = c(1,2))
image(.detangle_matrix(bsMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Jump Height", main = "Binary Segmentation")
image(.detangle_matrix(flMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Jump Height", main = "Fused Lasso")

####################################

bsMat <- form.matrix(bsFs_1JumpCI, resp.func = resp.length)
flMat <- form.matrix(flFs_1JumpCI, resp.func = resp.length)
zlim <- c(min(bsMat, flMat), max(bsMat, flMat))

#plot the unconditional length in absolute scale
par(mfrow = c(1,2))
image(.detangle_matrix(bsMat), zlim = zlim, xlab = "Jump Index", 
  ylab = "Jump Height", main = "Binary Segmentation")
contour(.detangle_matrix(bsMat), add = T, levels = 1, lwd = 3)
image(.detangle_matrix(flMat), zlim = zlim, xlab = "Jump Index", 
  ylab = "Jump Height", main = "Fused Lasso")
contour(.detangle_matrix(flMat), add = T, levels = 1, lwd = 3)

######################################

#plot the power (should be the same as pvalue)
bsMat <- form.matrix(bsFs_1JumpCI, resp.func = resp.power)
flMat <- form.matrix(flFs_1JumpCI, resp.func = resp.power)

par(mfrow = c(1,2))
image(.detangle_matrix(bsMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Jump Height", main = "Binary Segmentation")
image(.detangle_matrix(flMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Jump Height", main = "Fused Lasso")

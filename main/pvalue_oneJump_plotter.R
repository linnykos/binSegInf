rm(list=ls())
load("../results/pvalue_oneJump_2016-11-27.RData")

n <- 100
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
  mat  <- matrix (0, ncol = length(jump.loc), nrow = length (jump.height))
  colnames(mat) <- as.character(jump.loc)
  rownames(mat) <- as.character(jump.height)

  for (i in 1:length(res)){
    jump.percent  <-  as.numeric(strsplit(names(res)[i],split = "-")[[1]][3])
    true.loc  <- round(n*jump.percent)
  
    idx <- samp.selector(res[i], true.loc, ...)
    mat[i] <- resp.func(res[[i]], idx, true.loc)
  }
  
  mat
}

#some default resp.func
resp.power <- function(res, idx, true.loc){
  length(which(res[1,idx] <= alpha))/length(idx)
}
resp.exact <- function(res, idx, true.loc){
  length(idx)/ncol(res)
}
resp.variance <- function(res, idx, true.loc){
  mean(abs(res[2,idx]-true.loc))/n
}

bsMat <- form.matrix(bsFs_1JumpPValue, alpha)
flMat <- form.matrix(flFs_1JumpPValue, alpha)

#######################################

# unconditional plot with y-axis in absolute
par(mfrow = c(1,2))
image(.detangle_matrix(bsMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Jump Height", main = "Binary Segmentation")
contour(.detangle_matrix(bsMat), add = T, levels = 0.95, lwd = 3)

image(.detangle_matrix(flMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Jump Height", main = "Fused Lasso")
contour(.detangle_matrix(flMat), add = T, levels = 0.95, lwd = 3)

##############################

# unconditional plot with y-axis in log-scale
par(mfrow = c(1,2))
rownames(bsMat) <- as.character(log(as.numeric(rownames(bsMat))))
image(.detangle_matrix(bsMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "log Jump Height", main = "Binary Segmentation")
contour(.detangle_matrix(bsMat), add = T, levels = 0.95, lwd = 3)

rownames(flMat) <- as.character(log(as.numeric(rownames(flMat))))
image(.detangle_matrix(flMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "log Jump Height", main = "Fused Lasso")
contour(.detangle_matrix(flMat), add = T, levels = 0.95, lwd = 3)

####################################################

# percent chance of getting the right jump
bsMat <- form.matrix(bsFs_1JumpPValue, alpha, type = T, resp.func = resp.exact)
flMat <- form.matrix(flFs_1JumpPValue, alpha, type = T, resp.func = resp.exact)

par(mfrow = c(1,2))
image(.detangle_matrix(bsMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Jump Height", main = "Binary Segmentation")
contour(.detangle_matrix(bsMat), add = T, levels = 0.95, lwd = 3)

image(.detangle_matrix(flMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Jump Height", main = "Fused Lasso")
contour(.detangle_matrix(flMat), add = T, levels = 0.95, lwd = 3)

####################################################

#more lenient metric to getting the right jump
bsMat <- form.matrix(bsFs_1JumpPValue, alpha, resp.func = resp.variance)
flMat <- form.matrix(flFs_1JumpPValue, alpha, resp.func = resp.variance)
zlim = c(min(bsMat, flMat), max(bsMat, flMat))

rownames(bsMat) <- as.character(log(as.numeric(rownames(bsMat))))
rownames(flMat) <- as.character(log(as.numeric(rownames(flMat))))

par(mfrow = c(1,2))
image(.detangle_matrix(bsMat), zlim = zlim, xlab = "Jump Index", 
  ylab = "log Jump Height", main = "Binary Segmentation")
contour(.detangle_matrix(bsMat), add = T, levels = 100, lwd = 3)

image(.detangle_matrix(flMat), zlim = zlim, xlab = "Jump Index", 
  ylab = "log Jump Height", main = "Fused Lasso")
contour(.detangle_matrix(flMat), add = T, levels = 100, lwd = 3)


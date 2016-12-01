rm(list=ls())
load("../results/CI_oneJump_2016-11-30.RData")

n <- 100

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

form.matrix <- function(res, resp.func = resp.coverage, ...){
  mat  <- matrix (0, ncol = length(jump.loc), nrow = length (jump.height))
  colnames(mat) <- as.character(jump.loc)
  rownames(mat) <- as.character(jump.height)

  for (i in 1:length(res)){
    jump.percent  <-  as.numeric(strsplit(names(res)[i],split = "-")[[1]][3])
    true.loc  <- round(n*jump.percent)
  
    idx1 <- samp.selector(res[i], true.loc, ...)
    idx2 <- which(res[[i]][7,] == 0)
    idx <- intersect(idx1, idx2)
    mat[i] <- resp.func(res[[i]], idx, true.loc)
  }
  
  mat
}

resp.coverage <- function(res, idx, true.loc){
  length(intersect(which(res[2,idx] >= res[3,idx]), 
    which(res[2,idx] <= res[4,idx])))/length(idx)
}
resp.length <- function(res, idx, true.loc){
  len.vec <- apply(res[,idx], 2, function(x){x[4] - x[3]}); mean(len.vec)
}
resp.power <- function(res, idx, true.loc, n){
  length(unique(c(which(0 <= res[3,idx]), which(0 >= res[4,idx]))))/length(idx)
}


########################################

bsMat <- form.matrix(bsFs_1JumpCI)
flMat <- form.matrix(flFs_1JumpCI)

par(mfrow = c(1,2))
image(.detangle_matrix(bsMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Jump Height", main = "Binary Segmentation")
image(.detangle_matrix(flMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Jump Height", main = "Fused Lasso")


############################################

bsMat <- form.matrix(bsFs_1JumpCI, resp.func = resp.length)
flMat <- form.matrix(flFs_1JumpCI, resp.func = resp.length)
zlim <- c(min(bsMat, flMat), max(bsMat, flMat))

#plot the unconditional length in absolute scale
par(mfrow = c(1,2))
image(.detangle_matrix(bsMat), zlim = zlim, xlab = "Jump Index", 
  ylab = "Jump Height", main = "Binary Segmentation")
image(.detangle_matrix(flMat), zlim = zlim, xlab = "Jump Index", 
  ylab = "Jump Height", main = "Fused Lasso")

####################################################
#plot the unconditional length in log scale

rownames(bsMat) <- as.character(log(as.numeric(rownames(bsMat))))
rownames(flMat) <- as.character(log(as.numeric(rownames(flMat))))

par(mfrow = c(1,2))
image(.detangle_matrix(bsMat), zlim = zlim, xlab = "Jump Index", 
  ylab = "log Jump Height", main = "Binary Segmentation")
contour(.detangle_matrix(bsMat), add = T, levels = 3, lwd = 3)
contour(.detangle_matrix(bsMat), add = T, levels = 1.5, lwd = 3)

image(.detangle_matrix(flMat), zlim = zlim, xlab = "Jump Index", 
  ylab = "log Jump Height", main = "Fused Lasso")
contour(.detangle_matrix(flMat), add = T, levels = 3, lwd = 3)
contour(.detangle_matrix(flMat), add = T, levels = 1.5, lwd = 3)

#########################################33

#plot the power (should be the same as pvalue)
bsMat <- form.matrix(bsFs_1JumpCI, resp.func = resp.power)
flMat <- form.matrix(flFs_1JumpCI, resp.func = resp.power)

par(mfrow = c(1,2))
image(.detangle_matrix(bsMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Jump Height", main = "Binary Segmentation")
image(.detangle_matrix(flMat), zlim = c(0,1), xlab = "Jump Index", 
  ylab = "Jump Height", main = "Fused Lasso")

#######################################################

plot.intervals <- function(lis, limit = 50, ...){
  idx1 <- samp.selector(lis, true.loc, ...)
  res <- lis[[1]]
  idx2 <- which(res[7,] == 0)
  idx <- intersect(idx1, idx2)
  if(length(idx) > limit) idx <- idx[1:limit]
  
  plot(NA, xlim = c(1, length(idx)), ylim = c(min(res[3,idx]), max(res[4,idx])),
    xlab = "Index", ylab = "Value")
  for(i in 1:length(idx)){
    if(res[2,idx[i]] >= res[3,idx[i]] & res[2,idx[i]] <= res[4,idx[i]]){
      col = 4
    } else {
      col = 2
    }
    
    lines(x = rep(i,2), y = c(res[3,idx[i]], res[4,idx[i]]), col = col, lwd = 2)
    points(i, res[2,idx[i]], col = 1, cex = 2, pch = 16)
    points(i, res[1,idx[i]], col = col, cex = 1, pch = 16)
  }
  
  invisible()
}

## confidence intervals
par(mfrow = c(1,2))
plot.intervals(bsFs_1JumpCI[105])
title(main = "Binary Segmentation")

plot.intervals(flFs_1JumpCI[105])
title(main = "Fused Lasso")

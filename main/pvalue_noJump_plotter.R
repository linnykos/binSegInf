rm(list=ls())
load("../results/pvalue_noJump_2016-11-27.RData")

n <- 100
alpha <- 0.05
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

samp.selector <- function(lis, type = NA, func = function(x){x}){
  if(is.na(type)) return(1:ncol(lis[[1]]))

  bool.vec <- sapply(lis[[1]][2,], func)
  which(bool.vec)
}

idx <- samp.selector(bsFs_0JumpPValue[1])

##############################3

par(mfrow = c(1,2))
hist(bsFs_0JumpPValue[[1]][1,idx], freq = F, main = "Binary Segmentation", col = "gray", 
  breaks = 20, xlab = "p-value")
lines(density(bsFs_0JumpPValue[[1]][1,idx], cut = 0), col = "red", lwd = 2)

hist(flFs_0JumpPValue[[1]][1,idx], freq = F, main = "Fused Lasso", col = "gray", 
  breaks = 20, xlab = "p-value")
lines(density(flFs_0JumpPValue[[1]][1,idx], cut = 0), col = "red", lwd = 2)

############################

par(mfrow = c(1,2))
qqplot(bsFs_0JumpPValue[[1]][1,idx], seq(0, 1, length.out = ncol(bsFs_0JumpPValue[[1]][,idx])), 
  pch = 16, xlab = "Actual quantile", ylab = "Theoretical quantile")
lines(x = c(0,1), y = c(0, 1), col = "red", lwd = 2)

qqplot(flFs_0JumpPValue[[1]][1,idx], seq(0, 1, length.out = ncol(flFs_0JumpPValue[[1]][,idx])), 
  pch = 16, xlab = "Actual quantile", ylab = "Theoretical quantile")
lines(x = c(0,1), y = c(0, 1), col = "red", lwd = 2)

################################

graphics.off()
plot(bsFs_0JumpPValue[[1]][1,idx], flFs_0JumpPValue[[1]][1,idx], pch = 16, cex = 2,
  xlab = "Binary Segmentation", ylab = "Fused Lasso")

#########################################

mean(bsFs_0JumpPValue[[1]][3,])
mean(flFs_0JumpPValue[[1]][3,])

graphics.off()
par(mfrow = c(1,2))
plot(bsFs_0JumpPValue[[1]][3,], flFs_0JumpPValue[[1]][3,], pch = 16, cex = 2,
  xlab = "Binary Segmentation", ylab = "Fused Lasso", asp = T)
lines(c(0,10), c(0,10), col = "red", lwd = 2)

plot(sort(bsFs_0JumpPValue[[1]][3,]), pch = 16, cex = 2, xlab = "Index",
  ylab = "MSE", ylim = c(0, max(bsFs_0JumpPValue[[1]][3,])))
points(sort(flFs_0JumpPValue[[1]][3,]), pch = 16, cex = 2, col = 3)

####################################################

par(mfrow = c(1,3))
hist(bsFs_0JumpPValue[[1]][2,], col = "gray", breaks = 20)
hist(flFs_0JumpPValue[[1]][2,], col = "gray", breaks = 20)
plot(bsFs_0JumpPValue[[1]][2,], flFs_0JumpPValue[[1]][2,], pch = 16, cex = 2,
  xlab = "Binary Segmentation", ylab = "Fused Lasso")
lines(c(0,100), c(0,100), col = "red", lwd = 2)

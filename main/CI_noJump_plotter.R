rm(list=ls())
load("../results/CI_noJump_2016-11-27.RData")

n <- 100
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

samp.selector <- function(lis, type = NA, func = function(x){TRUE}){
  if(is.na(type)) return(1:ncol(lis[[1]]))

  bool.vec <- sapply(lis[[1]][2,], func)
  which(bool.vec)
}

compute.confidence <- function(mat){
  good.trials <- which(mat[6,] == 0)
  coverage.vec <- apply(mat[,good.trials], 2, function(x){
    if(0 >= x[2] & 0 <= x[3]) TRUE else FALSE
  })
  
  sum(coverage.vec)/ncol(mat[,good.trials])
}

compute.CIlength <- function(mat, keep.NA = F){
  if(keep.NA) good.trials <- 1:ncol(mat) else good.trials <- which(mat[6,] == 0)
  res <- apply(mat[,good.trials], 2, function(x){x[3]-x[2]})
  
  res[is.infinite(res)] <- NA
  res
}

############################################

compute.confidence(bsFs_0JumpCI[[1]])
compute.confidence(flFs_0JumpCI[[1]])

###########################################

bs.ci <- compute.CIlength(bsFs_0JumpCI[[1]])
fl.ci <- compute.CIlength(flFs_0JumpCI[[1]])
y.max <- max(bs.ci, fl.ci)

par(mfrow = c(1,2))
plot(sort(bs.ci), pch = 16, ylab = "CI width", ylim = c(0, y.max), main = "Binary Segmentation")
expected.diff <- qnorm(.975)-qnorm(.025)
lines(c(-1e5,1e5), rep(expected.diff, 2), col = "red", lwd = 2)

plot(sort(fl.ci), pch = 16, ylab = "CI width", ylim = c(0, y.max), main = "Fused Lasso")
expected.diff <- qnorm(.975)-qnorm(.025)
lines(c(-1e5,1e5), rep(expected.diff, 2), col = "red", lwd = 2)

############################################

graphics.off()
plot(sort(bs.ci), pch = 16, ylab = "CI width", ylim = c(0, y.max), main = "Binary Segmentation")
expected.diff <- qnorm(.975)-qnorm(.025)
lines(c(-1e5,1e5), rep(expected.diff, 2), col = "red", lwd = 2)

points(sort(fl.ci), pch = 16, col = 3)

#############################################

bs.ci2 <- compute.CIlength(bsFs_0JumpCI[[1]], T)
fl.ci2 <- compute.CIlength(flFs_0JumpCI[[1]], T)

idx <- sort(intersect(which(!is.na(bs.ci2)), which(!is.na(fl.ci2))))
sum(fl.ci2[idx] < bs.ci2[idx])/(length(fl.ci2[idx]))

par(mfrow = c(1,3))
plot(bs.ci2, fl.ci2, xlab = "Binary Segmentation CI Width",
  ylab = "Fused Lasso CI Width", pch = 16, cex = 2, asp = T)
lines(c(-1000,1000), c(-1000,1000), lwd = 2, col = "red")

hist(fl.ci2[idx] - bs.ci2[idx], col = "gray", breaks = 20, xlab = "FL - BS",
  main = "")

plot(sort(fl.ci2[idx] - bs.ci2[idx]), pch = 16, cex = 2,
  xlab = "Index", ylab = "FL - BS")

#####################################################

load("../results/pvalue_noJump_2016-11-27.RData")

par(mfrow = c(1,2))

plot(bs.ci2[1:500], bsFs_0JumpPValue[[1]][1,], pch = 16, cex = 2, main = "Binary Segmentation")
plot(fl.ci2[1:500], flFs_0JumpPValue[[1]][1,], pch = 16, cex = 2, main = "Fused Lasso")

#########################################################

plot.intervals <- function(lis, limit = 50, ...){
  idx1 <- samp.selector(lis, ...)
  res <- lis[[1]]
  idx2 <- which(res[6,] == 0)
  idx3 <- unique(which(0>=res[3,]), which(0<=res[2,]))
  idx <- setdiff(intersect(idx1, idx2), idx3)
  if(length(idx) > limit) idx <- idx[1:limit]
  
  plot(NA, xlim = c(1, length(idx)), ylim = c(min(res[2,idx]), max(res[3,idx])),
    xlab = "Index", ylab = "Value")
  for(i in 1:length(idx)){
    if(0 >= res[2,idx[i]] & 0 <= res[3,idx[i]]){
      col = 4
    } else {
      col = 2
    }
    
    lines(x = rep(i,2), y = c(res[2,idx[i]], res[3,idx[i]]), col = col, lwd = 2)
    points(i, 0, col = 1, cex = 2, pch = 16)
    points(i, res[1,idx[i]], col = col, cex = 1, pch = 16)
  }
  
  invisible()
}

par(mfrow = c(1,2))
plot.intervals(bsFs_0JumpCI)
title(main = "Binary Segmentation")

plot.intervals(flFs_0JumpCI)
title(main = "Fused Lasso")
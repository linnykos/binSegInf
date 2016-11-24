rm(list=ls())
load("../results/CI_noJump_bsFs_2016-11-18.RData")

res <- bsFs_0JumpCI
n <- 100
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

samp.selector <- function(lis, type = NA, func = function(x){TRUE}){
  if(is.na(type)) return(1:ncol(lis[[1]]))

  bool.vec <- sapply(lis[[1]][2,], func)
  which(bool.vec)
}

good.trials <- which(res[[1]][5,]==0)
idx <- intersect(samp.selector(res[1]), good.trials)

coverage.vec <- apply(res[[1]][,idx], 2, function(x){
  if(0 >= x[2] & 0 <= x[3]) TRUE else FALSE
})

sum(coverage.vec)/length(coverage.vec)

dist.diff <- apply(res[[1]][,idx], 2, function(x){x[3]-x[2]})
plot(sort(dist.diff), pch = 16, ylab = "CI width")
expected.diff <- qnorm(.975)-qnorm(.025)
lines(c(-1e5,1e5), rep(expected.diff, 2), col = "red", lwd = 2)

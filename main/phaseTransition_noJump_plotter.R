rm(list=ls())
load("../../Justin Kevin Max Ryan - binseginf/res/noJump-2016-09-11.RData")

par(mfrow = c(1,2))
hist(resFl[[1]], main = "Fused Lasso", col = rgb(.5,.5,.5))
hist(resBs[[1]], main = "Binary Seg", col = rgb(.5,.5,.5))
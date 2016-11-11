rm(list=ls())
load("../main/res/noJump-2016-09-11.RData")

res <- resFl1Jump
n <- 100
alpha <- 0.05
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

samp.selector <- function(lis, type = NA, func = function(x,y){x == y}){
  if(is.na(type)) return(1:ncol(lis[[1]]))
  
  jump.percent  <-  as.numeric(strsplit(names(lis), 
  split = "-")[[1]][3])
  true.loc  <- round(n*jump.percent)
  
  bool.vec <- sapply(lis[[1]][2,], func, y = true.loc)
  which(bool.vec)
}

idx <- samp.selector(res[1])
hist(res[[1]][1,idx], main = "P-values")
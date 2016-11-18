rm(list=ls())
load("../results/pvalue_oneJump_bsFs_2016-11-13.RData")

res <- bsFs_1JumpPValue
n <- 100
alpha <- 0.05

mat  <- matrix (0, ncol = length(jump.loc), nrow = length (jump.height))
colnames(mat) <- as.character(jump.loc)
rownames(mat) <- as.character(jump.height)

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

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
  two.sided <- sapply(res[[i]][1,], function(x){2*min(x, 1-x)})
  mat[i] <- length(which(two.sided <= alpha))/ncol(res[[i]][,idx])
}

image(.rotate(mat), zlim = c(0,1))
contour(.rotate(mat), add = T, levels = 0.95, lwd = 3)

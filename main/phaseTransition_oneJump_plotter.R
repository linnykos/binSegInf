rm(list=ls())
load("../../Justin Kevin Max Ryan - binseginf/res/oneJump-2016-09-15.RData")

n <- 200
MatFl  <- matrix (0, ncol = length(jump.loc), nrow = length (jump.height))
colnames(MatFl) <- as.character(jump.loc)
rownames(MatFl) <- as.character(jump.height)

for (i in 1:length (resFl1Jump)){
  jump.percent  <-  as.numeric(strsplit(names(resFl1Jump)[i], 
    split = "-")[[1]][3])
  true.loc  <- round(n*jump.percent)
  MatFl[i] <- length(which(resFl1Jump[[i]] == true.loc))/length(resFl1Jump[[i]])
}

image(t(MatFl)[,nrow(MatFl):1], zlim = c(0,1))

MatBs <- matrix (0, ncol = length(jump.loc), nrow = length (jump.height))
colnames(MatFl) <- as.character(jump.loc)
rownames(MatFl) <- as.character(jump.height)

for (i in 1:length (resFl1Jump)){
  jump.percent  <-  as.numeric(strsplit(names(resBs1Jump)[i], 
    split = "-")[[1]][3])
  true.loc  <- round(n*jump.percent)
  MatBs[i] <- length(which(resBs1Jump[[i]] == true.loc))/length(resBs1Jump[[i]])
}

image(t(MatBs)[,nrow(MatBs):1], zlim = c(0,1))


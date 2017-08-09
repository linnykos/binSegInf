rm(list=ls())
library(binSegInf)
set.seed(10)
## Simulation settings
n <- 20
numSteps <- 1
n_levs <- 1
lev <- 0
nsim <- 3000
sigma <- 1

## Generates one/two-jumped means
onejump <- function(lev,n){c(rep(0,n/2), rep(lev,n/2))}

## Main loop
pmat <- matrix(NA, ncol=2, nrow=nsim)
colnames(pmat) <- c("cp", "pv")

for(isim in 1:nsim){
  y <- onejump(lev,n) + rnorm(n,0,sigma)
  obj <- binSeg_fixedSteps(y, numSteps=numSteps)
  poly <- polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  pmat[isim,"cp"] <- obj$cp
  pmat[isim,"pv"] <- pvalue(y, poly, contrast)
  
  if(isim %% floor(nsim/10) == 0) cat('*')
}

plot(sort(pmat[,"pv"]), seq(0, 1, length.out = nsim))
lines(x = c(0,1), y = c(0,1), col = "red", lwd = 2)

res <- ks.test(pmat[,"pv"], punif, alternative = "two.sided")
res$p.value

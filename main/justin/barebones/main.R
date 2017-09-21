## Synopsis: produce noise-added randomized p-values for a simple algorithm that
## finds the single largest upward segment mean difference.
source("../main/justin/barebones/funs.R")
source("../main/justin/barebones/tests.R") ## unit tests

## Example settings
nsim = 500
sigma=1
sigma.add=.1
n=6
all.pvs4 = mclapply(1:nsim, function(isim){
    printprogress(isim,nsim)
    y = rep(0,n) + rnorm(n,0,sigma)
    added.noise = rep(0,n) + rnorm(n,0,sigma.add)
    rtg(y=y, sigma=sigma, shift=added.noise, sigma.add=sigma.add, nsim.inner=100)
}, mc.cores=4)
par(mfrow=c(5,1))
qqunif(unlist(all.pvs))
qqunif(unlist(c( all.pvs2)))
qqunif(unlist(c( all.pvs3)))
qqunif(unlist(c( all.pvs4)))
qqunif(unlist(c(all.pvs, all.pvs2,all.pvs3, all.pvs4)))

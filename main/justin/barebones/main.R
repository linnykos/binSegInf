## Synopsis: produce noise-added randomized p-values for a simple algorithm that
## finds the single largest upward segment mean difference.
source("../main/justin/barebones/funs.R")
source("../main/justin/barebones/tests.R") ## unit tests

## Example settings
nsim = 2000
sigma=1
sigma.add=.1
n=6
all.pvs = mclapply(1:nsim, function(isim){
    printprogress(isim,nsim)
    y = rep(0,n) + rnorm(n)
    rtg(y,sigma,sigma.add,nsim.inner=50)
}, mc.cores=3)
qqunif(unlist(all.pvs))

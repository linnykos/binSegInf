## Synopsis: produce noise-added randomized p-values for a simple algorithm that
## finds the single largest upward segment mean difference.
source("../main/justin/barebones/funs.R")
source("../main/justin/barebones/tests.R") ## unit tests

## Example settings
nsim = 10000
sigma=1
sigma.add=.1
n=6
all.pvs = mclapply(1:nsim, function(isim){
    printprogress(isim,nsim)
    y = rep(0,n) + rnorm(n,0,sigma)
    added.noise = rep(0,n) + rnorm(n,0,sigma.add)
    rtg(y=y, sigma=sigma, shift=added.noise, sigma.add=sigma.add, nsim.inner=100)
}, mc.cores=4)
qqunif(a..p.vs)



## Trying it over different signal sizes
nsim = 200
sigma = 1
sigma.add = 0.1
n = 6
all.pvs.list = list()
levs = c(0,1,3,5)
for(ii in 1:length(levs)){
    lev=levs[ii]
    cat("lev is", lev, fill=TRUE)
    all.pvs.list[[ii]] = mclapply(1:nsim, function(isim){
        printprogress(isim,nsim)
        y = c(rep(0,n/2), rep(lev,n/2)) + rnorm(n,0,sigma)
        added.noise = rep(0,n) + rnorm(n,0,sigma.add)
        rtg(y=y, sigma=sigma, shift=added.noise, sigma.add=sigma.add, nsim.inner=100)
    }, mc.cores=4)
}

## Plot it
cols = RColorBrewer::brewer.pal(length(levs),"Set1")
for(ii in 1:length(levs)){
    if(ii!=1) par(new=TRUE)
    qqunif(unlist(all.pvs.nonnull.list[[ii]]), col = cols[ii])
}

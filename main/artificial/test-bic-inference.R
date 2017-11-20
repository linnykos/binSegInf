
## Check uniformity
source("../main/artificial/artif-helpers.R")
n=10
mn = rep(0,n)
sigma=1
sigma.add = 0.2
nsim=100
results = lapply(1:nsim, function(isim){
    printprogress(isim,nsim)
    set.seed(isim)
    y = mn + rnorm(n,0,sigma)
    pvs = do_rfl_inference(y=y, max.numSteps=8,
                           consec=2, sigma=sigma, postprocess=TRUE,
                           locs=1:length(y), numIS=100, sigma.add = sigma.add, bits=1000)
    return(pvs)
    })
}, mc.cores=1)


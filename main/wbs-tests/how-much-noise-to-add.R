##' Synopsis: explore the effect of additive noise on the recovery and power of binseg


## Fixing a signal-to-noise ratio
n = 60
lev = 1
nsim = 1000
visc = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
numIntervals = round(seq(from=1/10,to=1.5,by=1/5)*n)
numSteps=4
locs = Map(function(my.numInterval)
    dosim_recovery(lev=lev,n=n,nsim=nsim,
                   numSteps=numSteps,
                   randomized=TRUE, numIS=100,
                   meanfun=fourjump,
                   numIntervals = my.numInterval,
                   mc.cores=mc.cores, locs=visc), numIntervals)

## Literally a prequel to the analysis with n=2000 example; timing things.
## n = length(y.orig)
library(microbenchmark)

datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))

## Generate nosie around the data

## Just do two types of inference: randwbs, randBS
ns = c(20,50,100,150,200, 250, 300)
rwbs.times = mclapply(1:length(ns), function(ii){
    ## Generate data
    n=ns[ii]
    printprogress(n, ns)
    sigma=1
    lev=2
    newmn = c(rep(0,n/2), rep(lev, n/2))
    return(microbenchmark({
        y = newmn + rnorm(n, 0, sigma)
        cumsum.y = cumsum(y)
        pvs.rwbs = do_rwbs_inference(y=y, max.numSteps=10, numIntervals=length(y),
                                     consec=2, sigma=sigma, postprocess=TRUE,
                                     better.segment=TRUE, locs=1:length(y),
                                     numIS=100, inference.type="pre-multiply",
                                     improve.nomass.problem=TRUE)
    }, times=10))
}, mc.cores = 7)


## Timing things
n=300
newmn = c(rep(0,n/2), rep(lev, n/2))
y = newmn + rnorm(n, 0, sigma)
Rprof("~/Desktop/timing-rwbs-inference.out")
pvs.rwbs = do_rwbs_inference(y=y, max.numSteps=4, numIntervals=length(y),
                             consec=2, sigma=sigma, postprocess=TRUE,
                             better.segment=TRUE, locs=1:length(y),
                             numIS=100, inference.type="pre-multiply",
                             improve.nomass.problem=TRUE)
Rprof(NULL)


pvs.rwbs = rep(NA,5)
microbenchmark({
    y = newmn + rnorm(n, 0, sigma)
    cumsum.y = cumsum(y)
source(file=file.path("../main/artificial/artif-helpers.R"))
    pvs.rwbs = do_rwbs_inference(y=y, max.numSteps=4, numIntervals=length(y),
                                 consec=2, sigma=sigma, postprocess=TRUE,
                                 better.segment=TRUE, locs=1:length(y),
                                 numIS=100, inference.type="pre-multiply",
                                 improve.nomass.problem=TRUE)
}, times=10))


getmedtime <- function(mytime){(mytime$time[3])/1000000000}
getmntime <- function(mytime){(mytime$time[2])/1000000000}
mntimes= sapply(rwbs.times, getmntime)
## mntimes[c(3,5,7)] = mntimes[c(3,5,78)]*1000
mntimes[c(3,5,7)] = mntimes[c(3,5,78)]*1000


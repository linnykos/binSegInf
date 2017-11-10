## Synopsis: one-jump inference examples
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")
outputdir = "../output"

## One-jump, one-step simulations
levs = c(0,1,2,3)
n=60
nsim=1000
numSteps=1

## Single simulation results
n = 60
numSteps = 10
nsim=1000
sigma=1
numIntervals=n



## Temporary run
levs = c(0,1,2,3)
n = 60
nsims=c(3000,700,500,250)

## nsims=500
## levs=0
n=60
numSteps=10
mc.cores=8
consec=2
visc = (n/2+((-1):1))
nsim=1
source("../main/wbs-tests/sim-helpers-with-stoptime.R")
nsim=200
nsims=c(5000,1000,500,250)
levs = c(0,1,2,3)
visc = (n/2+((-1):1))
## results = dosim_with_stoprule(lev=lev, n=n, nsim=nsim, numSteps=numSteps,randomized=TRUE,
##                               numIS=100, meanfun=onejump, mc.cores=4,
##                               inference.type="pre-multiply", consec=2, locs=visc)
results = Map(function(lev,nsim)dosim_with_stoprule(lev=lev,n=n,nsim=nsim,
                                                    numSteps=numSteps,
                                                    randomized=TRUE,numIS=100,
                                                    meanfun=onejump,
                                                    mc.cores=mc.cores,consec=consec,
                                                    locs=visc), levs, nsims)

## Save results
outputdir = "../output"
filename = "bic-wbs-onejump.Rdata"
save(list=c("results","levs","n","nsim","numSteps"), file=file.path(outputdir,filename))

## Synopsis: Check if onejump results are
source("../main/wbs-tests/sim-helpers.R")
levs = 0
n = 10
nsims=5000
numSteps=1
mc.cores=8
visc = (n/2+((-1):1))
## results = Map(function(lev,nsim) dosim(lev= lev,n=n,nsim=nsim,numSteps=numSteps,randomized=TRUE,numIS=100,
##                                        meanfun=onejump,mc.cores=mc.cores, locs=visc), levs, nsims)

lev=0
nsim=10000
results = dosim(lev= lev,n=n,nsim=nsim,numSteps=numSteps,randomized=TRUE,numIS=100,
      meanfun=onejump,mc.cores=mc.cores, locs=visc)

filename = "rand-wbs-onejump-one-step-nulls.Rdata"
save(list=c("results","levs","n","nsims","numSteps"), file=file.path(outputdir,filename))

## Synopsis: Try to make a better, improved segment test by using the winning
## intervals as segments

## Sample settings
set.seed(1)
sigma=1
lev=3
n=10
meanfun = onejump
numSteps=1
numIntervals=n
randomized=FALSE
visc = n/2 + c((-1),0,+1)
locs = visc

n=60
nsims=c(3000,700,500,250)
## numSteps=3
mc.cores=1
visc = (n/2+((-1):1))
## results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump, mc.cores=mc.cores)
source("../main/wbs-tests/sim-helpers.R")
results = Map(function(lev,nsim)dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps,
                                      randomized=randomized,numIS=100,meanfun=onejump,
                                      mc.cores=mc.cores, locs=visc, better.segment=TRUE), levs, nsims)



## On fourjump example
source("../main/wbs-tests/sim-helpers.R")
lev=1
n=60
mc.cores=4
meanfun = fourjump
numSteps=4
randomized=TRUE
visc = n/5 * (1:4)
mc.cores=4
nsim=1000
results = dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps, randomized=randomized,
                numIS=100, meanfun=meanfun, mc.cores=mc.cores, locs=visc,
                better.segment=TRUE, improve.nomass.problem = TRUE, min.num.things=20)
results.orig = dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps, randomized=randomize,numIS=100,
                     meanfun=meanfun, mc.cores=mc.cores, locs=visc, better.segment=FALSE,
                     improve.nomass.problem = TRUE, min.num.things=20)

filename = "better-segment-fourjump.Rdata"
save(list=c("results", "results.orig"), file=file.path(outputdir, filename))

## Load and plot better vs original segment test.
load(file=file.path(outputdir, filename))
qqunif(results$pv)
a = qqunif(results.orig$pv, plot.it=FALSE)
points(a,col='red')
legend("bottomright", legend=c("orig", "improved"), col = c("red", "black"), pch=c(16,16))





## ## On onejump example
## source("../main/wbs-tests/sim-helpers.R")
## outputdir = "../output"
## lev=1
## n=10
## mc.cores=4
## nsim=300
## meanfun = onejump
## numSteps=1
## randomized=TRUE
## ## visc = n/5 * (1:4)
## visc = n/2 + (c((-1):1))
## numIS=100
## improve.nomass.problem = TRUE
## results = dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps, randomized=randomized,numIS=numIS,
##                 meanfun=meanfun, mc.cores=mc.cores, locs=visc, better.segment=TRUE,
##                 improve.nomass.problem= improve.nomass.problem )
## results.orig = dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps, randomized=randomized,numIS=numIS,
##                      meanfun=meanfun, mc.cores=mc.cores, locs=visc, better.segment=FALSE,
##                      improve.nomass.problem =improve.nomass.problem )


## filename = "better-segment-onejump.Rdata"
## save(list=c("results.orig"), file=file.path(outputdir, filename))
## load(file=file.path(outputdir, filename))
## qqunif(results$pv)
## a = qqunif(results.orig$pv, plot.it=FALSE)
## points(a,col='blue')
## legend("bottomright", legend=c("orig", "improved"), col = c("blue", "black"), pch=c(16,16))


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
numSteps=1
mc.cores=8
visc = (n/2+((-1):1))
nsim=1
source("../main/wbs-tests/sim-helpers-with-stoptime.R")
nsim=40
pvs = dosim_with_stoprule(lev=0,n=10,nsim=nsim,numSteps= 5,randomized=TRUE,
                    numIS= 100,meanfun=onejump,mc.cores=1,
                    inference.type="pre-multiply", consec=2)

## Why is it still testing, when stoptime is zero? Is indexing off by onek?

## Why is the second value /always/ NaN?
results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps, randomized=FALSE, numIS=100, meanfun=onejump)


## Save results
outputdir = "../output"
filename = "fixed-wbs-onejump-one-step.Rdata"
save(list=c("results","levs","n","nsim","numSteps"), file=file.path(outputdir,filename))

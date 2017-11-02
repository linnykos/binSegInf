## Synopsis: one-jump inference examples
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")

## One-jump mean, one-step model
levs = c(0,1,2,3)
n=60
nsim=1000
numSteps=1
mc.cores=6
results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump, mc.cores=mc.cores)

nsim=1
numSteps=4
system.time({
dosim(lev=0, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump, mc.cores=1)
})

## Save results
outputdir = "../output"
filename = "rand-wbs-onejump-one-step.Rdata"
save(list=c("results","levs","n","nsim","numSteps"), file=file.path(outputdir,filename))


## One-jump mean, three-step model
levs = c(0,1,2,3)
n=60
nsim=1000
numSteps=3
mc.cores=6
results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump, mc.cores=mc.cores)

## Save results
outputdir = "../output"
filename = "rand-wbs-onejump-three-step.Rdata"
save(list=c("results","levs","n","nsim","numSteps"), file=file.path(outputdir,filename))



## One-jump mean, three-step model
levs = c(0,1,2,3)
n = 60
nsim = 1000
numSteps = 1
results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump, mc.cores=6)

## Save results
outputdir = "../output"
filename = "rand-wbs-onejump-three-step.Rdata"
save(list=c("results","levs","n","nsim","numSteps"), file=file.path(outputdir,filename))

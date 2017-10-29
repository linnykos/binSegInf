## Synopsis: one-jump inference examples
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")

## One-jump mean, one-step model
levs = c(0,1,2,3)
results = lapply(levs, dosim, n=60, nsim=1000, numSteps=1, randomized=TRUE, numIS=100, meanfun=onejump)

## Save results
outputdir = "../output"
filename = "fixed-wbs-onejump-one-step.Rdata"
save(list=c("results"), file=filename)


## One-jump mean, three-step model
levs = c(0,1,2,3)
results = lapply(levs, dosim, n=60, nsim=50000, numSteps=3, randomized=TRUE, numIS=100, meanfun=onejump)


## Save results
outputdir = "../output"
filename = "fixed-wbs-onejump-three-step.Rdata"
save(list=c("results"), file=filename)


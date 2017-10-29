## Synopsis: one-jump inference examples
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")

## One-jump simulations
levs = c(0,1,2,3)
results = lapply(levs, dosim, n=60, nsim=50000, numSteps=3, randomized=FALSE, numIS=100, meanfun=twojump)

## Save results
outputdir = "../output"
filename = "fixed-wbs-twojump-three-step.Rdata"
save(list=c("levs","pvs","locs","truths"), file=filename)

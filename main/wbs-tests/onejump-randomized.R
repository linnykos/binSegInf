## Synopsis: one-jump inference examples
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")

## One-jump simulations
levs = c(0,1,2,3)
onejumpsim(lev=0, n=60,nsim=1000,numSteps=3,randomized=TRUE, numIS=100)
results = lapply(levs, onejumpsim, n=60, nsim=50000, numSteps=3, randomized=FALSE, numIS=100)
pvs = lapply(results,function(a)a$pvs)
truths = results$truths

## Save results
outputdir = "../output"
filename = "fixed-wbs-onejump-three-step-visc.Rdata"
save(list=c("levs","pvs","locs","truths"), file=filename)

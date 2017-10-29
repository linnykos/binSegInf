## Synopsis: one-jump inference examples
source("../main/justin/plot-helpers.R")
source("../main/justin/sim-helpers.R")

## Three-jump simulations
levs = c(0,1,2,3)
results = lapply(levs, onejumpsim, n=60, nsim=50000, numSteps=1, randomized=TRUE)
pvs = lapply(results,function(a)a$pvs)
truths = results$truths

## Save results
outputdir = "../output"
filename = "fixed-wbs-onejump-one-step-visc.Rdata"
save(list=c("levs","pvs","locs","truths"), file=filename)



## Load results
load(file=file.path(outputdir,filename))

## Plot + Export results
pdf("nonrand-wbs.pdf")

## Process results
visc = (n/2+((-3):3))
visc = (n/2+((-1):1))
visc = (n/2)
locs = lapply(pvs, function(mypvs) as.numeric(names(mypvs)))
visc.pvs <- Map(function(mypv, myloc){
    mypv[myloc %in% (n/2+((-3):3))]},pvs, locs)
names(cond.pvs) = paste("jump-size=",levs)

qqunif(cond.pvs,cols=1:4)
maketitle(n=n, other="Three steps")

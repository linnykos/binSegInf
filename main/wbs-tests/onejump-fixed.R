## Synopsis: one-jump inference examples
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")
outputdir = "../output"

## One-jump, one-step simulations
levs = c(0,1,2,3)
n=60
nsim=10000
numSteps=1
visc = (n/2+((-1):1))
nsims=c(10000,4000,2000,1000)*10
mc.cores=4
results = Map(function(lev,nsim)dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps,randomized=FALSE,numIS=100,meanfun=onejump,mc.cores=mc.cores, locs=visc), levs, nsims)

## Save results
filename = "fixed-wbs-onejump-one-step.Rdata"
save(list=c("results","levs","n","nsim","numSteps"), file=file.path(outputdir,filename))


## One-jump, three-step simulations
levs = c(0,1,2,3)
n=60
numSteps=3
## results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps, randomized=FALSE, numIS=100, meanfun=onejump)
nsims=c(10000,4000,2000,1000)*10
visc = (n/2+((-1):1))
mc.cores=1
results = Map(function(lev, nsim)dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps,
                                       randomized=FALSE,numIS=100,
                                       meanfun=onejump,mc.cores=mc.cores,
                                       locs=visc), levs, nsims)

## Save results
filename = "fixed-wbs-onejump-three-step.Rdata"
save(list=c("results","levs","n","nsim","numSteps"), file=file.path(outputdir,filename))



## Load and extract results
filename = "fixed-wbs-onejump-one-step.Rdata"
load(file=file.path(outputdir,filename))
pvs = lapply(results, function(a)a[["pvs"]])
locs = lapply(pvs, function(mypvs) as.numeric(names(mypvs)))
visc = (n/2+((-1):1))


## Plot cond p-values over signal sizes
w = h = 5
pdf(file.path(outputdir,"fixed-wbs-onejump-one-step-visc.pdf"), width=w, height=h)
plot.visc(pvs,locs,visc,n,levs)
graphics.off()


## Load and extract results
filename = "fixed-wbs-onejump-three-step.Rdata"
load(file=file.path(outputdir,filename))
pvs = lapply(results, function(a)a[["pvs"]])
locs = lapply(results, function(a)a[["pvs"]])
truths = lapply(results,function(a)a$truth)
visc = (n/2+((-1):1))


## Plot null p-values
w = h = 5
pdf(file.path(outputdir,"fixed-wbs-onejump-three-step-nulls.pdf"), width=w, height=h)
plot.nulls(pvs,locs,truths,n,levs)
graphics.off()


## Plot cond visc p-values
pdf(file.path(outputdir,"fixed-wbs-onejump-three-step-visc.pdf"), width=w, height=h)
plot.visc(pvs,locs,visc,n,levs)
graphics.off()

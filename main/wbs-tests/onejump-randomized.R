## Synopsis: one-jump inference examples
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")

## One-jump mean, one-step model
levs = c(0,1,2,3)
n = 60
nsims=c(3000,700,500,250)
## nsims=500
## levs=0
numSteps=1
mc.cores=8
visc = (n/2+((-1):1))
results = Map(function(lev,nsim) dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps,randomized=TRUE,numIS=100,meanfun=onejump,mc.cores=mc.cores, locs=visc), levs, nsims)

## Save results
outputdir = "../output"
filename = "rand-wbs-onejump-one-step.Rdata"
save(list=c("results","levs","n","nsim","numSteps"), file=file.path(outputdir,filename))


## One-jump mean, three-step model
levs = c(0,1,2,3)
n=60
## nsim=1000
nsims=c(3000,700,500,250)
numSteps=3
mc.cores=6
visc = (n/2+((-1):1))
## results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump, mc.cores=mc.cores)
results = Map(function(lev,nsim)dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps,randomized=TRUE,numIS=100,meanfun=onejump,mc.cores=mc.cores, locs=visc), levs, nsims)

## Save results
outputdir = "../output"
filename = "rand-wbs-onejump-three-step.Rdata"
save(list=c("results","levs","n","nsim","numSteps"), file=file.path(outputdir,filename))





## Load and extract results
outputdir = "../output"
filename = "rand-wbs-onejump-three-step2.Rdata"
## filename = "rand-wbs-onejump-three-step.Rdata"
load(file=file.path(outputdir,filename))
pvs = lapply(results, function(a)a[["pvs"]])
locs = lapply(pvs, function(mypvs) as.numeric(names(mypvs)))
visc = (n/2+((-10):10))


## Plot cond p-values over signal sizes
w = h = 5
pdf(file.path(outputdir,"rand-wbs-onejump-one-step-visc.pdf"), width=w, height=h)
## pdf(file.path(outputdir,"rand-wbs-onejump-three-step-visc.pdf"), width=w, height=h)
plot.visc(pvs,locs,visc,n,levs)
graphics.off()


## Load and extract results
filename = "rand-wbs-onejump-three-step.Rdata"
load(file=file.path(outputdir,filename))
pvs = lapply(results, function(a)a[["pvs"]])
locs = lapply(results, function(a)a[["pvs"]])
truths = lapply(results,function(a)a$truth)
visc = (n/2+((-1):1))


## Plot null p-values
w = h = 5
pdf(file.path(outputdir,"rand-wbs-onejump-three-step-nulls.pdf"), width=w, height=h)
plot.nulls(pvs,locs,truths,n,levs)
graphics.off()


## Plot cond visc p-values
pdf(file.path(outputdir,"rand-wbs-onejump-three-step-visc.pdf"), width=w, height=h)
plot.visc(pvs,locs,visc,n,levs)
graphics.off()

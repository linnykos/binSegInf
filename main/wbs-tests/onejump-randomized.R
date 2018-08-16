## Synopsis: one-jump inference examples
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")
outputdir = "../output"

## One-jump mean, one-step model
levs = c(0,1,2,3)
n = 60
nsims=c(3000,700,500,250)*10
numSteps=1
mc.cores=8
visc = (n/2+((-1):1))
results = Map(function(lev,nsim) dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps,randomized=TRUE,numIS=100,meanfun=onejump,mc.cores=mc.cores, locs=visc), levs, nsims)

## Save results
filename = "rand-wbs-onejump-one-step.Rdata"
save(list=c("results","levs","n","nsims","numSteps"), file=file.path(outputdir,filename))



## One-jump mean, three-step model
## levs=0
## nsims=3000
n=60
levs = c(0,1,2,3)
nsims=c(3000,700,500,250)*3
numSteps=3
mc.cores=7
numIS = 100
visc = (n/2+((-1):1))
results = Map(function(lev,nsim)dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps,
                                      randomized=TRUE,numIS=numIS,meanfun=onejump,
                                      mc.cores=mc.cores, locs=visc, improve.nomass.problem=TRUE,
                                      inference.type="pre-multiply"), levs, nsims)

## Save results
filename = "rand-wbs-onejump-three-step-simpler.Rdata"
outputdir = "../output"
save(list=c("results","levs","n","nsims","numSteps"), file=file.path(outputdir,filename))


## Load and extract results
outputdir = "../output"
## outputdir = "~/Desktop"
## filename = "rand-wbs-onejump-one-step.Rdata"
## filename = "rand-wbs-onejump-three-step.Rdata"
load(file=file.path(outputdir,filename))
pvs = lapply(results, function(a)a[["pvs"]])
locs = lapply(pvs, function(mypvs) as.numeric(names(mypvs)))
visc = (n/2+((-10):10))


## Plot cond p-values over signal sizes
w = h = 5
pdf(file.path(outputdir,"rand-wbs-onejump-one-step-visc.pdf"), width=w, height=h)
## pdf(file.path(outputdir,"rand-wbs-onejump-three-step-visc.pdf"), width=w, height=h)
plot.visc(pvs,locs,visc,n,levs,mytitle="one jump + one step, near true jumps")
graphics.off()

pdf(file.path(outputdir,"rand-wbs-onejump-one-step-null.pdf"), width=w, height=h)
locs = lapply(pvs, function(mypvs) as.numeric(names(mypvs)))
truths = lapply(results,function(a)a$truth)
plot.nulls(pvs,locs,truths,n,levs, mytitle="one jump + one step, nulls")
graphics.off()



## Load and extract results
outputdir = "~/Desktop"
outputdir = "../output"
filename = "rand-wbs-onejump-three-step-moreIS.Rdata"
load(file=file.path(outputdir,filename))
pvs = lapply(results, function(a)a[["pvs"]])
locs = lapply(pvs, function(mypvs) as.numeric(names(mypvs)))
truths = lapply(results,function(a)a$truth)
visc = (n/2+((-1):1))


## Plot null p-values
w = h = 5
pdf(file.path(outputdir,"rand-wbs-onejump-three-step-nulls-higherIS.pdf"), width=w, height=h)
plot.nulls(pvs,locs,truths,n,levs, mytitle="one jump + three steps, nulls")
graphics.off()


## Plot cond visc p-values
pdf(file.path(outputdir,"rand-wbs-onejump-three-step-visc-higherIS.pdf"), width=w, height=h)
plot.visc(pvs,locs,visc,n,levs, mytitle="one jump + three steps, near true jumps")
graphics.off()

## Synopsis: one-jump inference examples
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")

## Two-jump, four-step model (RANDOMIZED)
levs = c(0,1,2,3)
n = 60
nsims=c(3000,700,500,250)
numSteps = 4
mc.cores = 6
visc = c((n/3+((-1):1)), (2*n/3+((-1):1)))
results = Map(function(lev,nsim) dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps,randomized=TRUE,numIS=100, meanfun=twojump,mc.cores=mc.cores, locs=visc), levs, nsims)

## Save results
outputdir = "../output"
filename = "rand-wbs-twojump-four-step.Rdata"
save(list=c("results","levs","n","nsim","numSteps"), file=file.path(outputdir,filename))




## Load and extract results
filename = "rand-wbs-twojump-four-step.Rdata"
load(file=file.path(outputdir,filename))
pvs = lapply(results, function(a)a[["pvs"]])
locs = lapply(pvs, function(mypvs) as.numeric(names(mypvs)))
visc = (n/2+((-1):1))
truths = lapply(results,function(a)a$truth)


## Plot cond p-values over signal sizes
w = h = 5
pdf(file.path(outputdir,"rand-wbs-twojump-four-step-visc.pdf"), width=w, height=h)
plot.visc(pvs,locs,visc,n,levs)
graphics.off()


## Plot null p-values
w = h = 5
pdf(file.path(outputdir,"rand-wbs-twojump-four-step-nulls.pdf"), width=w, height=h)
plot.nulls(pvs,locs,truths,n,levs)
graphics.off()

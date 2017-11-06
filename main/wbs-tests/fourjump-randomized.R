## Synopsis: one-jump inference examples
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")

## Four-jump, six-step simulations
levs = c(0,1,2,3)
n = 60
nsims=c(3000,700,500,250)
numSteps = 6
mc.cores = 8
visc = sapply(1:5, function(ii)ii*n/5+(-1):1)
results = Map(function(lev,nsim) dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps,randomized=TRUE,numIS=100,
                                       meanfun=fourjump,mc.cores=mc.cores, locs=visc), levs, nsims)

## Save results
outputdir = "../output"
filename = "rand-wbs-fourjump-six-step-visc.Rdata"
save(list=c("results","levs","n","nsim","numSteps"), file=file.path(outputdir,filename))


## Load and extract results
outputdir = "../output"
filename = "rand-wbs-fourjump-six-step.Rdata"
## filename = "rand-wbs-onejump-three-step.Rdata"
load(file=file.path(outputdir,filename))
pvs = lapply(results, function(a)a[["pvs"]])
locs = lapply(pvs, function(mypvs) as.numeric(names(mypvs)))
visc = (n/2+((-10):10))


## Plot cond p-values over signal sizes
w = h = 5
pdf(file.path(outputdir,"rand-wbs-fourjump-six-step-visc.pdf"), width=w, height=h)
## pdf(file.path(outputdir,"rand-wbs-onejump-three-step-visc.pdf"), width=w, height=h)
plot.visc(pvs,locs,visc,n,levs)
graphics.off()



## Synopsis: one-jump inference examples
## source("../main/wbs-tests/plot-helpers.R")
## source("../main/wbs-tests/sim-helpers.R")
outputdir = "../output"

## One-jump, one-step simulations
levs = c(0,1,2,3)
n=60
nsim=1000
numSteps = 10
n=60
numSteps=10
mc.cores=7
consec=2
visc = (n/2+((-1):1))
nsims=c(5000,1000,500,250)
levs = c(0,1,2,3)
visc = (n/2+((-1):1))
results = Map(function(lev,nsim)dosim_with_stoprule(lev=lev,n=n,nsim=nsim,
                                                    numSteps=numSteps,
                                                    randomized=TRUE,numIS=100,
                                                    meanfun=onejump,
                                                    mc.cores=mc.cores,consec=consec,
                                                    locs=visc), levs, nsims)

## Save results
outputdir = "../output"
filename = "bic-wbs-onejump.Rdata"
save(list=c("results","levs","n","nsims","visc", "numSteps"), file=file.path(outputdir,filename))


## Plot results
filename = "bic-wbs-onejump.Rdata"
load(file=file.path(outputdir,filename))
ilev = 1
results[[ilev]]$pvs
pvs = lapply(results, function(a)a[["pvs"]])
locs = lapply(pvs, function(mypvs) as.numeric(names(mypvs)))
visc = (n/2+((-1):1))


## Plot cond p-values over signal sizes
w = h = 5
pdf(file.path(outputdir,"wbs-onejump-stoprule-visc.pdf"), width=w, height=h)
plot.visc(pvs,locs,visc,n,levs)
graphics.off()



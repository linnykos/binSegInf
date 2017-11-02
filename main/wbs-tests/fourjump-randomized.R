## Synopsis: one-jump inference examples
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")

## Four-jump, six-step simulations
levs = c(0,1,2,3)
n = 60
nsim = 1000
numSteps = 6
mc.cores = 6
results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=fourjump, mc.cores=mc.cores)


## Save results
outputdir = "../output"
filename = "rand-wbs-fourjump-six-step-visc.Rdata"
save(list=c("results","levs","n","nsim","numSteps"), file=file.path(outputdir,filename))


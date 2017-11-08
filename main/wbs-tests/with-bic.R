## Synopsis: one-jump inference examples
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")
outputdir = "../output"

## One-jump, one-step simulations
levs = c(0,1,2,3)
n=60
nsim=1000
numSteps=1

## Single simulation results
n = 60
numSteps = 10
nsim=1000
sigma=1
numIntervals=n



## Temporary run
levs = c(0,1,2,3)
n = 60
nsims=c(3000,700,500,250)

## nsims=500
## levs=0
numSteps=1
mc.cores=8
visc = (n/2+((-1):1))
nsim=1
source("../main/wbs-tests/sim-helpers-with-stoptime.R")
nsim=20
results = dosim_with_stoprule(lev=0,n=10,nsim=nsim,numSteps= 10,randomized=TRUE,
                    numIS= 100,meanfun=onejump,mc.cores=4,
                    inference.type="pre-multiply", consec=2)


## Why is the second value /always/ NaN?
results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps, randomized=FALSE, numIS=100, meanfun=onejump)


## Save results
outputdir = "../output"
filename = "fixed-wbs-onejump-one-step.Rdata"
save(list=c("results","levs","n","nsim","numSteps"), file=file.path(outputdir,filename))






## Standalone code example

## mn = c(rep(0,n/2), rep(lev,n/2))
source("../main/wbs-tests/sim-helpers-with-stoptime.R")
nsim=20
lev = 0
n = 10
meanfun=onejump
sigma=1
numIntervals=n
mn = meanfun(lev,n)
set.seed(20)##21
y = mn + rnorm(n, 0, sigma)


## Fit initial WBS for a generous number of steps
numSteps=10
g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps,
                          inference.type='rows')

## Get ic-infused polyhedron
consec=2
ic_obj = get_ic(g$cp, g$y, consec=consec, sigma=sigma, type="bic")
ic_poly = ic_to_poly(ic_obj)
stoptime = ic_obj$stoptime

cp = g$cp[1:stoptime]
cp.sign = g$cp.sign[1:stoptime]
vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
v=vlist[[1]]
numIS=100
inference.type="pre-multiply"
randomize_wbsfs(v=v, winning.wbs.obj=g,
                sigma=sigma,
                numIS=numIS,
                inference.type=inference.type,
                cumsum.y=cumsum.y,
                cumsum.v=cumsum.v,
                stop.time=stoptime+consec,
                ic.poly=ic_poly)

nsim.is=200
pvs = sapply(vlist, function(v){
    cumsum.v = cumsum(v)
    return(randomize_wbsfs(v=v, winning.wbs.obj=g,
                           sigma=sigma,
                           numIS=numIS,
                           inference.type=inference.type,
                           cumsum.y=cumsum.y,
                           cumsum.v=cumsum.v,
                           stop.time=stoptime+consec,
                           ic.poly=ic_poly))
})



## Non-ic poly
source("../main/wbs-tests/sim-helpers-with-stoptime.R")
nsim=20
lev = 0
n = 10
meanfun=onejump
sigma=1
numIntervals=n
mn = meanfun(lev,n)
set.seed(1)##21
y = mn + rnorm(n, 0, sigma)

numSteps=1
g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps,
                          inference.type='rows')

nsim.is=200
stoptime = length(g$cp)
cp = g$cp[1:stoptime]
cp.sign = g$cp.sign[1:stoptime]
vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
pvs = sapply(vlist, function(v){
    cumsum.v = cumsum(v)
    return(randomize_wbsfs(v=v, winning.wbs.obj=g,
                           sigma=sigma,
                           numIS=numIS,
                           inference.type=inference.type,
                           cumsum.y=cumsum.y,
                           cumsum.v=cumsum.v))
})

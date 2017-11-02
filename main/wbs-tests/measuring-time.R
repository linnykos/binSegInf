library(microbenchmark )
## Load
n=60
lev=0
sigma = 1
## Generate some data
mn = onejump(lev,n)
set.seed(0)
y = mn + rnorm(n, 0, sigma)
numsteps=1

## Fit WBS
numIntervals = n
set.seed(1)
g.orig = wildBinSeg_fixedSteps(y, numIntervals=numIntervals,
                               numSteps= numSteps,
                               inference.type="rows")
poly = polyhedra(obj=g.orig$gamma, u=g.orig$u)
vlist <- make_all_segment_contrasts(g.orig)
numIS=100
v = vlist[[1]]
poly.pval2(y=y, poly=poly, v=v, sigma=sigma)$pv

cumsum.y = cumsum(y)
cumsum.v = cumsum(v)


## Rprof of inference
Rprof("~/Desktop/rand-wbs.out")
set.seed(0)
randomize_wbsfs(v, winning.wbs.obj=g.orig, numIS = 100, sigma=sigma,
                comprehensive=FALSE, inference.type="pre-multiply",
                cumsum.y=cumsum.y, cumsum.v=cumsum.v)
Rprof(NULL)
b1 = summaryRprof("~/Desktop/rand-wbs.out")

Rprof("~/Desktop/rand-wbs-orig.out")
set.seed(0)
randomize_wbsfs(v, winning.wbs.obj=g.orig, numIS = 100, sigma=sigma,
                comprehensive=FALSE, inference.type="rows",
                cumsum.y=cumsum.y, cumsum.v=cumsum.v)
Rprof(NULL)
b2= summaryRprof("~/Desktop/rand-wbs-orig.out")



## Microbenchmark of inference
c1 = microbenchmark({
    set.seed(0);
randomize_wbsfs(v, winning.wbs.obj=g.orig, numIS = 100, sigma=sigma,
                comprehensive=FALSE, inference.type="pre-multiply",
                cusum.y=cusum.y,cusum.v=cusum.v)}, times=5)

c2= microbenchmark({
    set.seed(0);
randomize_wbsfs(v, winning.wbs.obj=g.orig, numIS = 100, sigma=sigma,
                comprehensive=FALSE, inference.type="rows",
                cusum.y=cusum.y,cusum.v=cusum.v)}, times=5)


## Does the cusum_fast() speed up things more than cusum()?
n=60
lev=0
sigma = 1
## Generate some data
mn = onejump(lev,n)
set.seed(0)
y = mn + rnorm(n, 0, sigma)
cumsum.y = cumsum(y)
s=1
e=10
b=5
microbenchmark({cusum_fast(s=s,b=b,e=e,cumsums=cumsum.y)})
microbenchmark({cusum(s=s,b=b,e=e,y=y,contrast.vec=TRUE)})

## Is it actually the case that we are happy with the speed?

## Two-jump, four-step model
source("../main/wbs-tests/sim-helpers.R")
levs = c(0,1,2,3)
n = 60
nsim = 500
numSteps = 4
numIS=100
sigma=1
numIntervals=n
meanfun = twojump
numSteps=1
randomized=TRUE
nsim=1
set.seed(0)
numSteps=1
after.speedup1 = microbenchmark({
    dosim(lev=0, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE,
          numIS=100, meanfun=twojump, mc.cores=1, numIntervals=60,
          inference.type="pre-multiply")
}, times=3)

set.seed(0)
numSteps=3
after.speedup3 = microbenchmark({
    dosim(lev=0, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE,
          numIS=100, meanfun=twojump, mc.cores=1, numIntervals=60,
          inference.type="pre-multiply")
}, times=3)

set.seed(0)
numSteps=1
before.speedup1 = microbenchmark({
    dosim(lev=0, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100,
          meanfun=twojump, mc.cores=1, numIntervals=60, inference.type="rows")
},times=3)

set.seed(0)
numSteps=3
before.speedup3 = microbenchmark({
    dosim(lev=0, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100,
          meanfun=twojump, mc.cores=1, numIntervals=60, inference.type="rows")
},times=3)


Rprof("~/rand-wbs-before.out")
set.seed(0)
randomize_wbsfs(v, winning.wbs.obj=g.orig, numIS = 100, sigma=sigma,
                comprehensive=FALSE, inference.type="rows",
                cumsum.y=cumsum.y, cumsum.v=cumsum.v)
Rprof(NULL)
b1= summaryRprof("~/rand-wbs-before.out")
head(b1$by.total, 20)
head(b1$by.self, 20)

## Rprof of inference
Rprof("~/rand-wbs-after.out")
set.seed(0)
randomize_wbsfs(v, winning.wbs.obj=g.orig, numIS = 100, sigma=sigma,
                comprehensive=FALSE, inference.type="pre-multiply",
                cumsum.y=cumsum.y, cumsum.v=cumsum.v)
Rprof(NULL)
b2 = summaryRprof("~/rand-wbs-after.out")
head(b2$by.total, 20)
head(b2$by.self, 20)


head(b1$by.total, 20)
head(b2$by.total, 20)


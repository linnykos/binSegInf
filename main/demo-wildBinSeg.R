## Fixed threshold wild binary segmentation

## Generate data
n=20
lev=5
numIntervals = 10
thresh = 1
nsim.is=10
set.seed(0)
cpv = CpVector(n=n, jump.height=c(0,lev),jump.loc=1/2)
y = cpv$data


#######################
## WBS-fixed-thresh ###
#######################

## Fix method and pfun
method = wildBinSeg_fixedThresh
pfun = poly.pval
pfun2 = randomized_wildBinSeg_pv

## Run method /once/, collect things
obj = method(y,thresh, numIntervals=numIntervals, verbose=FALSE)

## Collect a polyhedron
poly <- polyhedra(obj)

## Make contrasts
vs = make_all_segment_contrasts(obj)

## Do inference and return
pv = sapply(vs, function(v){
    pfun(y=y,v=v, G=poly$gamma, u=poly$u, sigma=sigma)$pv})

## Optionally, compare it directly to a different pfun2
pv2 = sapply(vs, function(v){ pfun2(y=y,v=v,poly=poly,sigma=sigma,
                                    nsim.is=nsim.is, numIntervals=numIntervals)})

#######################
## WBS-fixed-steps ###
#######################

## Test settings
n = 10
lev = 5
mu = c(rep(0,n/2),rep(lev,n/2))
sigma=1
set.seed(0)
y <- mu + rnorm(n)
numSteps = 3
set.seed(0)
numIntervals=10
intervals = generate_intervals(length(y), numIntervals)

## Run method /once/, collect things
numSteps=5
obj = wildBinSeg_fixedSteps(y, numSteps, intervals=intervals, verbose=TRUE)

## Fix method and pfun
pfun = poly.pval
pfun2 = randomized_wildBinSeg_pv


## Apply IC rule, get stoptime and polyhedra for it.
ic_wrapper(obj, y=y, sigma=sigma, returntype="polyhedra")
ic_wrapper(obj, y=y, sigma=sigma, returntype="stoptime")


## Collect a polyhedron
poly <- polyhedra(obj)

## Make contrasts
vs = make_all_segment_contrasts(obj)

## Do inference and return
pv = sapply(vs, function(v){
    pfun(y=y,v=v, G=poly$gamma, u=poly$u, sigma=sigma)$pv})

## Optionally, compare it directly to a different pfun2
pv2 = sapply(vs, function(v){ pfun2(y=y,v=v,poly=poly,sigma=sigma,
                                    nsim.is=nsim.is, numIntervals=numIntervals)})

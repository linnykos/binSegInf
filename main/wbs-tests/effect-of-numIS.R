lev = levs=0
n=10
## nsims=c(3000,700,500,250)
nsims=1000
numSteps=3
mc.cores=8
numIS = 200
visc = (n/2+((-1):1))
## results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump, mc.cores=mc.cores)
mc.cores=1
## if(isim==63)browser()
isim=63
set.seed(isim)
printprogress(isim, nsim)
sigma=1
numIntervals=n

## Generate some data
source("../main/wbs-tests/sim-helpers.R")
meanfun=onejump
mn = meanfun(lev,n)
y = mn + rnorm(n, 0, sigma)
cumsum.y=cumsum(y)

## Fit WBS
locs = visc = (n/2+((-1):1))
g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
poly = polyhedra(obj=g$gamma, u=g$u)
vlist <- make_all_segment_contrasts_from_cp(cp=g$cp, cp.sign=g$cp.sign, n=n)

## Retain only the guys we want
retain = which((g$cp %in% locs))
if(length(retain)==0){
    return(list(pvs=c(), null.true=c()))
}

## Get the p-values
vlist = vlist[retain] ## Added
v = vlist[[1]]
cumsum.v=cumsum(v)
nsim=20
inference.type = "pre-multiply"


## numIS.list = c(100,200,500,1000,2000)
numIS.list = c(3000,5000,100000)
nsim2 = 1000

pvslist = lapply(numIS.list, function(numIS){
    print(numIS)
    pvs = mclapply(1:nsim2, function(isim){
        printprogress(isim,nsim2)
        randomize_wbsfs(v=v, winning.wbs.obj=g, sigma=sigma,
                        numIS=numIS, inference.type=inference.type,
                        cumsum.y=cumsum.y,cumsum.v=cumsum.v)
    }, mc.cores=6)
    print("")
    pvs=unlist(pvs)
    return(pvs)
})

pvslist.sup2000 = pvslist
## save(list="pvslist.sub2000", file=file.path(outputdir, "effect-of-numIS-sub2000.Rdata"))
save(list="pvslistl.sup2000", file=file.path(outputdir, "effect-of-numIS-sup2000.Rdata"))


outputdir = "~/Desktop"
filename="effect-of-numIS.pdf"
pdf(file=file.path(outputdir, filename), width=15,height=5)
par(mfrow=c(1,5))
for(jj in 1:5){
    hist(pvslist[[jj]],main="", xlab="")
    title(main=paste("numIS=", numIS.list[jj]))
    title(sub=paste("n=10, numSteps=3, lev=0"))
}
graphics.off()

##' Synopsis: see the distribution of p-values (within a single choice of y and
##' v) as we vary numIS
lev = levs=0
n=10
## nsims=c(3000,700,500,250)
nsims=1000
numSteps=3
mc.cores=8
numIS = 200
visc = (n/2+((-1):1))
## results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump, mc.cores=mc.cores)
mc.cores=1
## if(isim==63)browser()
isim=63
set.seed(isim)
printprogress(isim, nsim)
numIntervals=n

## Generate some data
source("../main/wbs-tests/sim-helpers.R")
meanfun=onejump
mn = meanfun(lev,n)
y = mn + rnorm(n, 0, sigma)
cumsum.y=cumsum(y)

## Fit WBS
locs = visc = (n/2+((-1):1))
g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
poly = polyhedra(obj=g$gamma, u=g$u)
vlist <- make_all_segment_contrasts_from_cp(cp=g$cp, cp.sign=g$cp.sign, n=n)

## Retain only the guys we want
retain = which((g$cp %in% locs))
if(length(retain)==0){
    return(list(pvs=c(), null.true=c()))
}

## Get the p-values
vlist = vlist[retain] ## Added
v = vlist[[1]]
cumsum.v=cumsum(v)
nsim=20
inference.type = "pre-multiply"

## numIS.list = c(100,200,500,1000,2000)
numIS.list = c(3000,4000)
nsim2 = 500
pvslist = list()
for(ii in 1:2){
    pvslist[[ii]] = mclapply(1:nsim2, function(isim){
        printprogress(isim,nsim2)
        randomize_wbsfs(v=v, winning.wbs.obj=g, sigma=sigma,
                        numIS=numIS.list[ii], inference.type=inference.type,
                        cumsum.y=cumsum.y,cumsum.v=cumsum.v)
    }, mc.cores=3)
}
save(pvslist, file.path(outputdir,"effect-of-numIS-sub2000.Rdata"))0wb




## Incorporating the improve-nomass-problem; do I get uniform null p-values now?
lev = levs=0
n=10
## nsims=c(3000,700,500,250)
nsims=1000
numSteps=3
mc.cores=7
numIS = 200
visc = (n/2+((-1):1))
isim=63
sigma=1
numIntervals=n

## Generate some data
source("../main/wbs-tests/sim-helpers.R")
meanfun=onejump
nsim=1000
results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps,
                 randomized=TRUE, numIS=100, meanfun=onejump, mc.cores=mc.cores,
                 improve.nomass.problem=improve.nomass.problem)


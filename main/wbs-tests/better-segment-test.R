## Synopsis: Try to make a better, improved segment test by using the winning
## intervals as segments

## Sample settings
set.seed(1)
sigma=1
lev=3
n=10
meanfun = onejump
numSteps=1
numIntervals=n
randomized=FALSE
visc = n/2 + c((-1),0,+1)
locs = visc

n=60
nsims=c(3000,700,500,250)
## numSteps=3
mc.cores=1
visc = (n/2+((-1):1))
## results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump, mc.cores=mc.cores)
source("../main/wbs-tests/sim-helpers.R")
results = Map(function(lev,nsim)dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps,
                                      randomized=randomized,numIS=100,meanfun=onejump,
                                      mc.cores=mc.cores, locs=visc, better.segment=TRUE), levs, nsims)



## On fourjump example
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
## whichlevs=1:2
## whichlevs = 3:5
levs = c(0,1,2,3,4)[whichlevs]
n = 200
meanfun = fourjump
numSteps = 4
randomized = TRUE
## visc = n/5
## visc = 3*n/5 + c((-2):2)
visc = 2*n/5 + c((-2):2)
mc.cores = 8
min.num.things = 10
bits=1000
## nsims = seq.int(from=500, to=100, length=5)[whichlevs]
nsims = c(1500, 750, 300, 100, 100)[whichlevs]
source("../main/wbs-tests/sim-helpers.R")
for(ilev in 1:length(whichlevs)){
    lev = levs[ilev]
    nsim = nsims[ilev]
    results = dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps, randomized=randomized,
                    numIS=100, meanfun=meanfun, mc.cores=mc.cores, locs=visc,
                    better.segment=TRUE, improve.nomass.problem = TRUE, inference.type="pre-multiply",
                    min.num.things=min.num.things, bits=bits)
    results.orig = dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps,
                         randomized=randomized,numIS=100, meanfun=meanfun,
                         mc.cores=mc.cores, locs=visc, better.segment=FALSE,
                         inference.type = "pre-multiply",
                         improve.nomass.problem = TRUE, min.num.things=min.num.things,
                         bits=bits)

    ## Save the result
    if(lev == 0.5) lev= "onehalf" ## For saving purposes
    ## filename= paste0("better-segment-fourjump-lev", lev, ".Rdata")
    filename= paste0("better-segment-fourjump-80-lev", lev, ".Rdata")
    save(list=c("results", "results.orig"), file=file.path(outputdir, filename))
    print(paste0("saved to ", filename))
}

## Load and make QQ plot better vs original segment test.
lev=3
filename= paste0("better-segment-fourjump-lev", lev, ".Rdata")
load(file=file.path(outputdir, filename))
qqunif(results$pv)
a = qqunif(results.orig$pv, plot.it=FALSE)
points(a,col='red')
legend("bottomright", legend=c("orig", "improved"), col = c("red", "black"), pch=c(16,16))

## Make power table
levs = c(0,2,3,4)
powermat = matrix(ncol=2,nrow=5)
colnames(powermat) = c("power", "powers.orig")
rownames(powermat) = levs
for(lev in levs){
    print(lev)
    ## Load and clean
    filename = paste0("better-segment-fourjump-lev", lev, ".Rdata")
    load(file=file.path(outputdir, filename))
    pvs.orig = results.orig$pv
    pvs.orig = pvs.orig[!is.na(pvs.orig)]
    pvs = results$pvs
    pvs = pvs[!is.na(pvs)]

    ## Calculate original powers
    power.orig = sum(pvs.orig < 0.05/4)/length(pvs.orig)
    power = sum(pvs< 0.05/4)/length(pvs)
    powermat[lev+1,] = c(power, power.orig)
}


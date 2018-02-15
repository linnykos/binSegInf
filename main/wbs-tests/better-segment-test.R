## Synopsis: Try to make a better, improved segment test by using the winning
## intervals as segments

## On fourjump example
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
whichlevs = 4:6
## whichlevs = 1:3
levs = c(0,0.5,1,2,3,4)[whichlevs]
n = 200
meanfun = fourjump
numSteps = 4
randomized = TRUE
visc = 3*n/5 + c((-2):2)
## visc = 2*n/5 + c((-2):2)
print(visc)
mc.cores = 8
min.num.things = 10
bits=1000
nsims = c(1500, 1000, 750, 300, 100, 100)[whichlevs]
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
    filename= paste0("better-segment-fourjump-120-lev", lev, ".Rdata")
    ## filename= paste0("better-segment-fourjump-80-lev", lev, ".Rdata")
    save(list=c("results", "results.orig"), file=file.path(outputdir, filename))
    print(paste0("saved to ", filename))
}

## Load and make QQ plot better vs original segment test.
## lev = 3
## filename= paste0("better-segment-fourjump-lev", lev, ".Rdata")
## load(file=file.path(outputdir, filename))
## qqunif(results$pv)
## a = qqunif(results.orig$pv, plot.it=FALSE)
## points(a,col='red')
## legend("bottomright", legend=c("orig", "improved"), col = c("red", "black"), pch=c(16,16))

## Make power table
levs = c(0,0.5,1,2,3,4)
powermat = matrix(ncol=2,nrow=6)
colnames(powermat) = c("power", "powers.orig")
rownames(powermat) = levs
for(ilev in 1:6){
    lev = levs[ilev]
    if(lev == 0.5) lev= "onehalf" ## For saving purposes
    ## Load and clean
    filename = paste0("better-segment-fourjump-80-lev", lev, ".Rdata")
    ## filename = paste0("better-segment-fourjump-lev", lev, ".Rdata")
    load(file=file.path(outputdir, filename))
    pvs.orig = results.orig$pv
    pvs.orig = pvs.orig[!is.na(pvs.orig)]
    pvs = results$pvs
    pvs = pvs[!is.na(pvs)]

    ## Calculate original powers
    power.orig = sum(pvs.orig < 0.05/4)/length(pvs.orig)
    power = sum(pvs< 0.05/4)/length(pvs)
    powermat[ilev,] = c(power, power.orig)
}
round(powermat,3)
xtable::xtable(powermat,digits=3)




## Why is this thing so better at lev1? Load the two up, and compare seed by seed.
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
source("../main/wbs-tests/sim-helpers.R")
lev=1
nsim=750
filename = paste0("better-segment-fourjump-lev", lev, ".Rdata")
load(file=file.path(outputdir, filename))
pvs.orig = results.orig$pv
pvs.orig = pvs.orig[!is.na(pvs.orig)]

## The first one is already huely
pvs = results$pvs
pvs = pvs[!is.na(pvs)]

plot(pvs~pvs.orig)


source("../main/wbs-tests/sim-helpers.R")
nsim=1
randomized=TRUE
mc.cores=1
meanfun=fourjump
numSteps=4
n=200
visc = 3*n/5 + c((-2):2)
bits=1000
min.num.things = 10
nsim = 5
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

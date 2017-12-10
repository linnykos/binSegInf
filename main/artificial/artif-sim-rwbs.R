# Data directory
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))

##' Multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1, n){
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}

## Simulation settings
onesim_rwbs <- function(y.orig){

    ## Add bootstrapped residuals around a cleaned mean, with known sigma
    y = newmn + bootstrap_sample(resid.cleanmn)

    ## The idea is to add bootstrap noise, then see the 1.
    ## conditional/unconditional powers of the tests done at the viscinities of
    ## the true guys, and 2. uniformity (or lack thereof) in the tests conducted
    ## regarding null locations.

    par(mfrow=c(3,3))
    plot(y)
    lines(cleanmn,col='red')
    for(ii in 1:8){
        plot(cleanmn + bootstrap_sample(y-cleanmn))
        lines(cleanmn,col='red')
    }

    ## Just do two types of inference: randwbs, randBS
    pvs.rwbs = do_rwbs_inference(y=y, max.numSteps=10, numIntervals=length(y),
                                 consec=2, sigma=sigma, postprocess=TRUE,
                                 better.segment=TRUE, locs=1:length(y),
                                 numIS=100, inference.type="pre-multiply",
                                 improve.nomass.problem=TRUE, bits=1000)
    return(pvs.rwbs)
}
nsim = 50
results = mclapply(1:nsim,function(isim){
    return(onesim_rwbs())
},mc.cores=8)

## Reading and seeing speed from file
read.time.from.file <- function(myfile){
    sort(readLines(myfile))
}


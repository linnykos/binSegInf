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
onesim_rbs <- function(y.orig, bits=1000){

    ## Add bootstrapped residuals around a cleaned mean, with known sigma
    y = newmn + bootstrap_sample(resid.cleanmn)

    ## The idea is to add bootstrap noise, then see the 1.
    ## conditional/unconditional powers of the tests done at the viscinities of
    ## the true guys, and 2. uniformity (or lack thereof) in the tests conducted
    ## regarding null locations.

    ## Just do two types of inference: randwbs, randBS
    pvs.rbs = do_rbs_inference(y=y, max.numSteps=10, numIntervals=length(y),
                               consec=2, sigma=sigma, postprocess=TRUE,
                               better.segment=TRUE, locs=1:length(y), numIS=100,
                               inference.type="pre-multiply",
                               improve.nomass.problem=TRUE, bits=bits,
                               write.file=TRUE)
    return(pvs.rbs)
}

nsim = 50
bits = 5000
results = mclapply(1:nsim, function(isim){
    printprogress(isim,nsim)
    return(onesim_rbs(y.orig, bits=bits))
},mc.cores=8)

## Reading and seeing speed from file
read.time.from.file <- function(myfile){
    sort(readLines(myfile))
}

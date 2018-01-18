# Data directory
datadir = "../data"
outputdir = "../output"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))

## Simulation driver
onesim_rwbs <- function(y.orig){

    ## Add bootstrapped residuals around a cleaned mean, with known sigma
    ## estimated in-sample.
    sigma = sd(y.orig[1:200])
    y = newmn[-(1:200)] + bootstrap_sample(resid.cleanmn[-(1:200)])

    ## The idea is to add bootstrap noise, then see the 1.
    ## conditional/unconditional powers of the tests done at the viscinities of
    ## the true guys, and 2. uniformity (or lack thereof) in the tests conducted
    ## regarding null locations.

    ## Just do two types of inference: randwbs, randBS
    pvs.rwbs = do_rwbs_inference(y=y, max.numSteps=10, numIntervals=length(y),
                                 consec=2, sigma=sigma, postprocess=TRUE,
                                 better.segment=TRUE, locs=1:length(y),
                                 numIS=100, inference.type="pre-multiply",
                                 improve.nomass.problem=TRUE, bits=1000, write.time=TRUE)
    return(pvs.rwbs)
}

nsim = 300 # 50
results = mclapply(1:nsim,function(isim){
    printprogress(isim,nsim)
    return(onesim_rwbs(y.orig))
},mc.cores=5)


save(results, file=file.path(outputdir, "artif-rwbs.Rdata"))

## Reading and seeing speed from file
read.time.from.file <- function(myfile){
    sort(readLines(myfile))
}


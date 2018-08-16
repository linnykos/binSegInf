## Synopsis: Want to see if marginalized WBS has p-value equal to one

outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
dosim_compare_rwbs_to_wbs <- function(n=20, numIntervals=n, nsim, lev=1,
                                      mc.cores=4, type = c("wbs.nonrand", "wbs.rand"),
                                      visc=1:n, meanfun, bits=1000,
                                      min.num.things=30, verbose=FALSE,seed=TRUE){

    ## visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
    type = match.arg(type)
    start.time = Sys.time()
    pvs.list = mclapply(1:nsim, function(isim){
        if(!is.null(seed)) set.seed(isim)
        printprogress(isim, nsim,
                      lapsetime = round(difftime(Sys.time(), start.time,
                                            units = "secs"), 2))
        pvs = dosim_compare(type=type, n=n, lev=lev, numIntervals=numIntervals,
                            numIS=100, meanfun=meanfun, visc=visc, numSteps=4,
                            max.numIS=2000, bits=bits,
                            min.num.things, verbose=verbose)
        if(any(pvs==1)) stop(msg=paste0("seed is ", isim))
        return(pvs)
    }, mc.cores=mc.cores)
    return(pvs.list)
}

## Isolating a cases in which the pvalue is equal to 1.
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")

## Setting for running this code on a server
meanfun = fourjump_spiky
n = 200
numIntervals = n
mc.cores=8
nsims = seq(from=3000,to=1500,length=length(levs))/nchunk
visc.firsthalf = 1:(n/2)
bits=2000
min.num.things = 30

## Actual simulation code for the server (just change ichunk=1,2,3)
lev = 0.25
meanfun = fourjump_spiky
nsim=2000
## Marginalized WBS
myresult = dosim_compare_rwbs_to_wbs(numIntervals=numIntervals, n=n, lev=lev,
                                     nsim=nsim, mc.cores=mc.cores, type = "wbs.rand",
                                     meanfun=meanfun,
                                     visc=visc.firsthalf,
                                     bits=bits,
                                     min.num.things=min.num.things)

## In that case, what is actually going on? How can we fix it? (See things)
isim=129123
pvs = dosim_compare(type="wbs.rand", n=n, lev=lev, numIntervals=numIntervals,
                    numIS=100, meanfun=meanfun, visc=visc, numSteps=4,
                    max.numIS=2000, bits=bits,
                    min.num.things, verbose=verbose)


## How long does it take for the first non-one case to occur? (What factors
## matter? bits?)
isim=129123
pvs = dosim_compare(type="wbs.rand", n=n, lev=lev, numIntervals=numIntervals,
                    numIS=100, meanfun=meanfun, visc=visc, numSteps=4,
                    max.numIS=200000000, bits=bits,
                    min.num.things, verbose=verbose)

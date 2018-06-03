3## Synopsis: Want to see if marginalization of WBS inference is worth it, in
## terms of power! Rerunning experiments for different-shaped signal
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


## Setting for running this code on a server
meanfun = fourjump_spiky
n = 200
nchunk = 3
numIntervals = n
mc.cores = 4
levs = c(0.5, 1)
nsims = seq(from=3000,to=1500,length=length(levs))/nchunk
visc.firsthalf = 1:(n/2)
bits=2000
min.num.things = 30

## Actual simulation code for the server (just change ichunk=1,2,3)
ichunk = 3
results.wbs = results.rwbs = list()
for(ilev in 1:length(levs)){
    lev = levs[ilev]
    printprogress(lev, levs, fill=TRUE)
    nsim = nsims[ilev]

    ## Spiky
    filename = paste0("wbs-vs-rwbs-spiky-chunk-", ichunk, ".Rdata")
    meanfun = fourjump_spiky

    ## ## Hybrid
    ## filename = paste0("wbs-vs-rwbs-hybrid-chunk-", ichunk, ".Rdata")
    ## meanfun = fourjump_hybrid

    ## ## Same size
    ## filename = paste0("wbs-vs-rwbs-samesize-chunk-", ichunk, ".Rdata")
    ## meanfun = fourjump_samesize

    ## Marginalized WBS
    cat(fill=TRUE)
    cat("Marginalized", fill=TRUE)
    results.rwbs[[ilev]] = dosim_compare_rwbs_to_wbs(numIntervals=numIntervals, n=n, lev=lev,
                                                     nsim=nsim, mc.cores=mc.cores, type = "wbs.rand",
                                                     meanfun=meanfun,
                                                     visc=visc.firsthalf,
                                                     bits=bits,
                                                     min.num.things=min.num.things)
    save(results.rwbs, results.wbs, file=file.path(outputdir, filename))

    cat("Nonmarginalized", fill=TRUE)
    ## Non-Marginalized WBS
    results.wbs[[ilev]] = dosim_compare_rwbs_to_wbs(numIntervals=numIntervals, n=n, lev=lev,
                                                    nsim=nsim, mc.cores=mc.cores, type = "wbs.nonrand",
                                                    meanfun=meanfun,
                                                    visc=visc.firsthalf,
                                                    bits=bits,
                                                    min.num.things=min.num.things)
    save(results.rwbs, results.wbs, file=file.path(outputdir, filename))
}

## Aggregate things across chunks (code this up when done.)
extract.from.results <- function(ilev, type){
    nchunk=3
    all.things = lapply(1:nchunk,function(ichunk){
        filename = paste0("wbs-vs-rwbs-spiky-chunk-", ichunk, ".Rdata")
        load(file=file.path(outputdir, filename))
        if(type=='rwbs'){
            results.rwbs[[ilev]]
        } else {
            results.wbs[[ilev]]
        }
    })
    do.call(c, all.things)
}

rwbs.lev1.results = extract.from.results(1, "rwbs")
rwbs.lev2.results = extract.from.results(2, "rwbs")
wbs.lev1.results = extract.from.results(1, "wbs")
wbs.lev2.results = extract.from.results(2, "wbs")


lev1.pvs.list = list(
    rwbs.lev1.pvs = (unlist(sapply(rwbs.lev1.results, function(a)a$pvs))),
    wbs.lev1.pvs = (unlist(sapply(wbs.lev1.results, function(a)a$pvs)))
)

lev2.pvs.list = list(
    rwbs.lev2.pvs = (unlist(sapply(rwbs.lev2.results, function(a)a$pvs))),
    wbs.lev2.pvs = (unlist(sapply(wbs.lev2.results, function(a)a$pvs)))
)

## Get rid of 1's, temporarily
lev1.pvs.list[[1]] = lev1.pvs.list[[1]][which(lev1.pvs.list[[1]]!=1)]
lev2.pvs.list[[1]] = lev2.pvs.list[[1]][which(lev2.pvs.list[[1]]!=1)]

## Make qqplots
filename = "fourjump-spike-qqplot-lev1.jpg"
jpeg(file=file.path(outputdir, filename), width=500,height=500)
qqunif(lev1.pvs.list, cols=c(1,2))
graphics.off()


filename = "fourjump-spike-qqplot-lev2.jpg"
jpeg(file=file.path(outputdir, filename), width=500,height=500)
qqunif(lev2.pvs.list, cols=c(1,2))
graphics.off()



## Nonuniform p-values, but
## allresults.rwbs.lev1 = list()
## allresults.rwbs[[ilev]] = c(allresults.rwbs[[ilev]], results.rwbs[[ilev]])

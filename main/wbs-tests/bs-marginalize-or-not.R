## Synopsis: See if we should marginalize out the additive noise in binary
## segmentation. i.e. look at the power. NBS means Noise-added BS. RBS means
## randomized (marginalized) additive noise BS.

outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
dosim_compare_rbs_to_nbs<- function(n=20, numIntervals=n, nsim, lev=1,
                                    mc.cores=4, type = c("sbs.rand.nonmarg", "sbs.rand", "sbs.nonrand"),
                                    visc=1:n, meanfun, bits=1000,
                                    min.num.things=30, verbose=FALSE,
                                    seed=NULL, sigma.add=.2, sigma=1){

    type = match.arg(type)
    start.time = Sys.time()
    pvs.list = mclapply(1:nsim, function(isim){
        if(!is.null(seed)) set.seed(isim)
        printprogress(isim, nsim,
                      lapsetime = round(difftime(Sys.time(), start.time,
                                            units = "secs"), 2))
        pvs = dosim_compare(type=type, n=n, lev=lev, numIntervals=numIntervals,
                            numIS=10, meanfun=meanfun, visc=visc, numSteps=4,
                            max.numIS=2000, bits=bits,
                            min.num.things=min.num.things, verbose=verbose,
                            sigma.add=sigma.add, sigma=sigma)
        ## if(any(pvs==1)) stop(msg=paste0("seed is ", isim))
        return(pvs)
    }, mc.cores=mc.cores)
    return(pvs.list)
}

## Setting for running this code on a server
n = 20
nchunk = 3
numIntervals = n
mc.cores = 7
levs = c(0,3)
nsims = c(1000,1000)
bits=5000
min.num.things = 20
sigma=1
meanfun = fourjump

## Actual simulation code for the server (just change ichunk=1,2,3)
ichunk = 3
filename = paste0("rbs-vs-nbs-chunk-", ichunk, ".Rdata")
results.rbs = results.nbs = list()
for(ilev in 1:length(levs)){
    lev = levs[ilev]
    printprogress(lev, levs, "levels", fill=TRUE)
    nsim = nsims[ilev]

    cat(fill=TRUE)
    cat("Nonmarginalized", fill=TRUE)
    results.rbs[[ilev]] = dosim_compare_rbs_to_nbs(numIntervals=numIntervals,
                                                   n=n, lev=lev, nsim=nsim,
                                                   mc.cores=mc.cores,
                                                   type = "sbs.rand.nonmarg",
                                                   meanfun=meanfun,
                                                   bits=bits,
                                                   min.num.things=min.num.things,
                                                   sigma.add=sigma*0.2)
    save(results.rbs, results.nbs, file=file.path(outputdir, filename))
}

for(ilev in 1:length(levs)){
    lev = levs[ilev]
    printprogress(lev, levs, "levels", fill=TRUE)
    nsim = nsims[ilev]

    cat(fill=TRUE)
    cat("Marginalized", fill=TRUE)
    results.nbs[[ilev]] = dosim_compare_rbs_to_nbs(numIntervals=numIntervals,
                                                   n=n, lev=lev, nsim=nsim,
                                                   mc.cores=mc.cores,
                                                   type = "sbs.rand",
                                                   meanfun=meanfun,
                                                   bits=bits,
                                                   min.num.things=min.num.things,
                                                   sigma.add=sigma*0.2)
    save(results.rbs, results.nbs, file=file.path(outputdir, filename))
}


## Trying plain
mc.cores=4
nsim=2000
    lev= 3
    plain = dosim_compare_rbs_to_nbs(numIntervals=numIntervals,
                                                   n=n, lev=lev, nsim=nsim,
                                                   mc.cores=mc.cores,
                                                   type = "sbs.rand",
                                                   meanfun=meanfun,
                                                   bits=bits,
                                                   min.num.things=min.num.things,
                                                   sigma.add=0)
    ## save(plain, file=file.path(outputdir, "rbs-vs-nbs-plain.Rdata"))


## Load and make qqplots.
outputdir="../output"
ichunk=1
filename = paste0("rbs-vs-nbs-chunk-", ichunk, ".Rdata")
load(file=file.path(outputdir, filename))
pvs.nbs.flat.master = lapply(1:3, function(ichunk){
    filename = paste0("rbs-vs-nbs-chunk-", ichunk, ".Rdata")
    load(file=file.path(outputdir, filename))
    sapply(results.nbs[[1]], function(a)a[,"pvs"])
})
pvs.rbs.flat.master = lapply(1:3, function(ichunk){
    filename = paste0("rbs-vs-nbs-chunk-", ichunk, ".Rdata")
    load(file=file.path(outputdir, filename))
    sapply(results.rbs[[1]], function(a)a[,"pvs"])
})
pvs.nbs.nonflat.master = lapply(1:3, function(ichunk){
    filename = paste0("rbs-vs-nbs-chunk-", ichunk, ".Rdata")
    load(file=file.path(outputdir, filename))
    sapply(results.nbs[[2]], function(a)a[,"pvs"])
})
pvs.rbs.nonflat.master = lapply(1:3, function(ichunk){
    filename = paste0("rbs-vs-nbs-chunk-", ichunk, ".Rdata")
    load(file=file.path(outputdir, filename))
    sapply(results.rbs[[2]], function(a)a[,"pvs"])
})

filename = "bs-marginalize-or-not.jpg"
jpeg(file=file.path(outputdir, filename), width=1000, height=500)

## Overall p-value distribution
par(mfrow=c(1,2))
qqunif(list(marginalized=unlist(pvs.nbs.nonflat.master),
            nonmarginalized=unlist(pvs.rbs.nonflat.master),
            plain=unlist(lapply(plain,function(a)a$pvs))),
       cols=c("black", "red","blue"))
title(main="All p-values")




## Locations
locs.nbs.nonflat.master = lapply(1:3, function(ichunk){
    filename = paste0("rbs-vs-nbs-chunk-", ichunk, ".Rdata")
    load(file=file.path(outputdir, filename))
    sapply(results.nbs[[1]], function(a)a[,"locs"])
})
locs.rbs.nonflat.master = lapply(1:3, function(ichunk){
    filename = paste0("rbs-vs-nbs-chunk-", ichunk, ".Rdata")
    load(file=file.path(outputdir, filename))
    sapply(results.rbs[[1]], function(a)a[,"locs"])
})

## Only p-values at correct locations
locs = c(4,8,12,16)
nbs.pvs.correct = unlist(pvs.nbs.nonflat.master)[which(unlist(locs.nbs.nonflat.master) %in% locs)]
rbs.pvs.correct = unlist(pvs.rbs.nonflat.master)[which(unlist(locs.rbs.nonflat.master) %in% locs)]
locs.plain = unlist(lapply(plain, function(a)a$locs))
pvs.plain = unlist(lapply(plain, function(a)a$pvs))
plain.correct = pvs.plain[which(locs.plain %in% locs)]
qqunif(list(marginalized = nbs.pvs.correct,
            nonmarginalized = rbs.pvs.correct,
            plain = plain.correct), cols = c("black", "red","blue"))
title(main="p-values of tests at correct locations")
graphics.off()


## Data plot
dataplotfilename
jpeg(file=file.path(outputdir, dataplotfilename), width=500,height=500)


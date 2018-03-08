## Synopsis: Newly writing the script for comparing things, from the bottom up.
## Some issues addressed first (erase these as you go).

library(genlassoinf)
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
onecompare <- function(lev=0, nsim=1000, mc.cores=8, meanfun=onejump, visc=NULL, numSteps=1, bits=50, n=60, numIS=200,
                       max.numIS=1000){

    all.names = c("fl.rand", "fl.rand.plus", "fl.nonrand", "sbs.rand",
                  "sbs.nonrand", "wbs.rand", "wbs.nonrand", "cbs.rand",
                  "cbs.nonrand")

    all.results = lapply(all.names, function(myname){
        print(myname)
        mclapply(1:nsim, function(isim) {
            printprogress(isim,nsim);
            dosim_compare(type=myname, n=n, lev=lev, numIS=numIS, meanfun=meanfun, visc=visc, numSteps=numSteps, bits=bits,
                          max.numIS=max.numIS)
        }, mc.cores=mc.cores)
    })
    names(all.results) = all.names
    return(all.results)
}
## Synopsis: Compare all methods

library(genlassoinf)
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
onecompare <- function(lev=0, nsim=1000, mc.cores=8, meanfun=onejump, visc=NULL, numSteps=1, bits=50, n=60, numIS=200,
                       max.numIS=1000, verbose=FALSE){

    all.names = c("fl.rand", "fl.rand.plus", "fl.nonrand", "sbs.rand",
                  "sbs.nonrand", "wbs.rand", "wbs.nonrand", "cbs.rand",
                  "cbs.nonrand")

    all.results = lapply(all.names, function(myname){
        print(myname)
        mclapply(1:nsim, function(isim) {
            printprogress(isim,nsim);
            dosim_compare(type=myname, n=n, lev=lev, numIS=numIS, meanfun=meanfun, visc=visc, numSteps=numSteps, bits=bits,
                          max.numIS=max.numIS, verbose=verbose)
        }, mc.cores=mc.cores)
    })
    names(all.results) = all.names
    return(all.results)
}

jj = 1
whichlev.list = list(1, 2, 3, 4, 5:6, 7:8, 9:10, 11:12)
whichlev = whichlev.list[[jj]]
levs = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4)
mc.cores = 8
whichlev =
nsims = c(3000,3000,3000,3000, seq(from=3000,to=1000,length=5),
            round(seq(from=600, to=300, length=4) ))[whichlev]
n = 200
visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
print(levs)
results.by.lev = list()
for(ilev in 1:length(levs)){
    printprogress(mylev, levs, "levels running")
    cat(fill=TRUE)
    results.by.lev[[ whichlev[ilev] ]] = onecompare(lev=levs[ilev],
                                                    nsim=nsims[ilev], meanfun=fourjump, visc=visc.fourjump,
                                                    numSteps=4, bits=3000, mc.cores=mc.cores, n=n, numIS=100,
                                                    max.numIS=2000, verbose=verbose)
    filename = paste0("compare-methods-fourjump-final-", ilev, ".Rdata")
    save(list=c("results.by.lev","levs","nsim", "n"), file=file.path(outputdir, filename))
    print(filename)
}

nsim=1
ilev=1
myname = "cbs.rand"
source("../main/wbs-tests/sim-helpers.R")
a = dosim_compare(type=myname, n=n, lev=lev, numIS=numIS, meanfun=meanfun,
                  visc=1:n, numSteps=numSteps, bits=bits, max.numIS=max.numIS, verbose=TRUE)


mc.cores = 4
nsims = 4
levs = 1
ilev=1
visc.fourjump = 1:n
verbose=FALSE
a =  onecompare(lev=levs[ilev], nsim=nsims[ilev], meanfun=fourjump,
                visc=visc.fourjump, numSteps=4, bits=3000, mc.cores=mc.cores,
                n=n, numIS=100, max.numIS=100, verbose=verbose)
## Check that all of them work
myname = "sbs.rand"
bits=3000
max.numIS = 100
source("../main/wbs-tests/sim-helpers.R")
a = dosim_compare(type=myname, n=n, lev=lev, numIS=numIS, meanfun=meanfun,
                  visc=1:n, numSteps=numSteps, bits=bits, max.numIS=max.numIS,
                  sigma.add)


## Before I start running this, what is stopping me

## Synopsis: See randomized tests power compared to the nonrandomized version
library(genlassoinf)
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
onecompare <- function(lev=0, filename, nsim=1000, mc.cores=8, meanfun=onejump, visc=NULL, numSteps=1){

    result.fl.rand = result.fl.nonrand = result.sbs.rand =
    result.sbs.nonrand = result.wbs.rand = result.wbs.nonrand =
    result.cbs.rand = result.cbs.nonrand = NA

    print("fl.rand")
    result.fl.rand =  mclapply(1:nsim, function(isim) {
        printprogress(isim,nsim);
        dosim_compare(type="fl.rand", n=60, lev=lev, numIS=200, meanfun=meanfun, visc=visc, numSteps=numSteps)
    }, mc.cores=mc.cores)

    print("fl.nonrand")
    result.fl.nonrand =  mclapply(1:nsim, function(isim) {
        printprogress(isim,nsim);
        dosim_compare(type="fl.nonrand", n=60, lev=lev, numIS=200, meanfun=meanfun, visc=visc, numSteps=numSteps)
    }, mc.cores=mc.cores)

    print("sbs.rand")
    result.sbs.rand =  mclapply(1:nsim, function(isim) {
        printprogress(isim,nsim);
        dosim_compare(type="sbs.rand", n=60, lev=lev, numIS=200, meanfun=meanfun, visc=visc, numSteps=numSteps)
    }, mc.cores=mc.cores)

    print("sbs.nonrand")
    result.sbs.nonrand =  mclapply(1:nsim, function(isim) {
        printprogress(isim,nsim);
        dosim_compare(type="sbs.nonrand", n=60, lev=lev, numIS=200, meanfun=meanfun, visc=visc, numSteps=numSteps)
    }, mc.cores=mc.cores)

    print("wbs.nonrand")
    result.wbs.nonrand =  mclapply(1:nsim, function(isim) {
        printprogress(isim,nsim);
        dosim_compare(type="wbs.nonrand", n=60, lev=lev, numIS=200, meanfun=meanfun, visc=visc, numSteps=numSteps)
    }, mc.cores=mc.cores)

    print("wbs.rand")
    result.wbs.rand =  mclapply(1:nsim, function(isim) {
        printprogress(isim,nsim);
        dosim_compare(type="wbs.rand", n=60, lev=lev, numIS=200, meanfun=meanfun, visc=visc, numSteps=numSteps)
    }, mc.cores=mc.cores)

    print("cbs.rand")
    result.cbs.rand =  mclapply(1:nsim, function(isim) {
        printprogress(isim,nsim);
        dosim_compare(type="cbs.rand", n=60, lev=lev, numIS=200, meanfun=meanfun, visc=visc, numSteps=numSteps)
    }, mc.cores=mc.cores)

    print("cbs.nonrand")
    result.cbs.nonrand =  mclapply(1:nsim, function(isim) {
        printprogress(isim,nsim);
        dosim_compare(type="cbs.nonrand", n=60, lev=lev, numIS=200, meanfun=meanfun, visc=visc, numSteps=numSteps)
    }, mc.cores=mc.cores)

    save(list=c("result.fl.rand", "result.fl.nonrand", "result.sbs.rand",
                "result.sbs.nonrand", "result.wbs.rand", "result.wbs.nonrand",
                "result.cbs.rand", "result.cbs.nonrand"),
         file=filename)
}

## onecompare(lev=1, filename=file.path(outputdir,'compare-methods-lev1-onejump-with-cbs-halve-fix.Rdata'),
##            nsim=10000, meanfun=onejump, visc=c(28:32), numSteps=1)
visc.fourjump = unlist(lapply(c(1:4)*60/5, function(ii) ii+c(-1,0,1)))
onecompare(lev=0, filename=file.path(outputdir,'compare-methods-lev0-fourjump.Rdata'), nsim=5000, meanfun=fourjump, visc=visc.fourjump, numSteps=4)
onecompare(lev=1, filename=file.path(outputdir,'compare-methods-lev1-fourjump.Rdata'), nsim=5000, meanfun=fourjump, visc=visc.fourjump, numSteps=4)
onecompare(lev=2, filename=file.path(outputdir,'compare-methods-lev2-fourjump.Rdata'), nsim=5000, meanfun=fourjump, visc=visc.fourjump, numSteps=4)
onecompare(lev=3, filename=file.path(outputdir,'compare-methods-lev3-fourjump.Rdata'), nsim=5000, meanfun=fourjump, visc=visc.fourjump, numSteps=4)
graphics.off()



## Run the nonnull simulation
onecompare(lev=3, filename=file.path(outputdir,'compare-methods-lev3.Rdata'))

## Load and plot
## load(file=file.path(outputdir, 'compare-methods-lev0.Rdata'))
load(file=file.path(outputdir, 'compare-methods-lev1.Rdata'))
load(file=file.path(outputdir, 'compare-methods-lev1-fourjump.Rdata'))
## load(file=file.path(outputdir, 'compare-methods-lev2.Rdata')) #
## load(file=file.path(outputdir, 'compare-methods-lev3.Rdata'))
## pdf(file=file.path(outputdir,"compare-methods-lev0.pdf"), width=5, height=5)
## pdf(file=file.path(outputdir,"compare-methods-lev1.pdf"), width=5, height=5)
pdf(file=file.path(outputdir,"compare-methods-lev1-fourjump.pdf"), width=5, height=5)
## pdf(file=file.path(outputdir,"compare-methods-lev2.pdf"), width=5, height=5)
## pdf(file=file.path(outputdir,"compare-methods-lev3.pdf"), width=5, height=5)

all.results = list(result.fl.nonrand, result.fl.rand, result.sbs.nonrand,
                   result.sbs.rand, result.wbs.nonrand, result.wbs.rand,
                   result.cbs.nonrand, result.cbs.rand)[c(3,4,5,6,7,8)]

all.names = c("fl.nonrand", "fl.rand", "sbs.nonrand",
                "sbs.rand", "wbs.nonrand", "wbs.rand",
                "cbs.nonrand", "cbs.rand")[c(3,4,5,6,7,8)]

pvs.list = Map( function(myresult){
    myresult = lapply(myresult, function(a){colnames(a) = c("pvs", "locs");a})
    pvs = do.call(rbind, myresult)[,1]
    ## tab = tab[which(apply(tab,1,function(myrow)!all(is.na(myrow)))),]
    pvs = pvs[!is.na(pvs)]
    return(pvs)
}, all.results)
pvs.list = lapply(pvs.list, unlist)
names(pvs.list) = all.names
mar = c(4.5,4.5,2.5,0.5)
## cols = RColorBrewer::brewer.pal(8,"Set1")
cols = rep(RColorBrewer::brewer.pal(3,"Set2")[c(1,2,3)],each=2)
qqunif_line(pvs.list, cols=cols, names=all.names, lty=c(1,2,1,2,1,2),lwds=rep(4,6))
title(main=expression(delta==0))
graphics.off()


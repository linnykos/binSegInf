## Synopsis: See randomized tests power compared to the nonrandomized version
library(genlassoinf)
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
onecompare <- function(lev=0, nsim=1000, mc.cores=8, meanfun=onejump, visc=NULL, numSteps=1, bits=50, n=60, numIS=200){

    all.names = c("fl.rand", "fl.nonrand", "sbs.rand",
                 "sbs.nonrand", "wbs.rand", "wbs.nonrand",
                 "cbs.rand", "cbs.nonrand")

    all.results = lapply(all.names, function(myname){
        print(myname)
        mclapply(1:nsim, function(isim) {
            printprogress(isim,nsim);
            dosim_compare(type=myname, n=n, lev=lev, numIS=numIS, meanfun=meanfun, visc=visc, numSteps=numSteps, bits=bits)
        }, mc.cores=mc.cores)
    })
    names(all.results) = all.names
    return(all.results)
}

levs = c(0, 0.5, 1, 1.5, 2)[4:5]
results.by.lev = list()
mc.cores=8
nsim=3000
nsims=seq(from=3000,to=1000,length=5)[4:5]
n=200
visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
for(ilev in 1:length(levs)){
    mylev = levs[ilev]
    nsim = nsims[ilev]
    print(mylev)
    results.by.lev[[ilev]] = onecompare(lev=mylev,
                                        nsim=nsim, meanfun=fourjump, visc=visc.fourjump,
                                        numSteps=4, bits=1000, mc.cores=mc.cores, n=200, numIS=200)
    ## save(list=c("results.by.lev","levs","nsim", "n"), file="compare-methods-fourjump-45.Rdata")
    save(list=c("results.by.lev","levs","nsim", "n"), file="compare-methods-fourjump-123.Rdata")
}








## load(file=file.path(outputdir, 'compare-methods-lev3.Rdata'))
## pdf(file=file.path(outputdir,"compare-methods-lev3.pdf"), width=5, height=5)

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


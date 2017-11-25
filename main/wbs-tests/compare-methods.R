## Synopsis: See randomized tests power compared to the nonrandomized version
library(genlassoinf)
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
onecompare <- function(lev=0, nsim=1000, mc.cores=8, meanfun=onejump, visc=NULL, numSteps=1, bits=50, n=60, numIS=200){

    all.names = c("fl.rand", "fl.nonrand", "sbs.rand",
                 "sbs.nonrand", "wbs.rand", "wbs.nonrand",
                 "cbs.rand", "cbs.nonrand")[1:2]

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

## levs = c(0, 0.5, 1, 1.5, 2)[4:5]
## levs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)[6:7]
## levs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)[8:9]
levs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
results.by.lev = list()
mc.cores=8
nsim=3000
nsims=c(seq(from=3000,to=1000,length=5), round(seq(from=600, to=300, length=4) ))[8:9]
n=200
visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
print(levs)
for(ilev in 1:length(levs)){
    mylev = levs[ilev]
    nsim = nsims[ilev]
    print(mylev)
    results.by.lev[[ilev]] = onecompare(lev=mylev,
                                        nsim=nsim, meanfun=fourjump, visc=visc.fourjump,
                                        numSteps=4, bits=1000, mc.cores=mc.cores, n=200, numIS=200)
    ## save(list=c("results.by.lev","levs","nsim", "n"), file=file.path(outputdir,"compare-methods-fourjump-45.Rdata"))
    ## save(list=c("results.by.lev","levs","nsim", "n"), file=file.path(outputdir, "compare-methods-fourjump-123.Rdata"))
    ## save(list=c("results.by.lev","levs","nsim", "n"), file=file.path(outputdir, "compare-methods-fourjump-67.Rdata"))
    ## save(list=c("results.by.lev","levs","nsim", "n"), file=file.path(outputdir, "compare-methods-fourjump-89.Rdata"))
    save(list=c("results.by.lev","levs","nsim", "n"), file=file.path(outputdir, "compare-methods-fourjump-only-fl.Rdata"))
}



## Load 123
outputdir = "../output"

## Aggregate the results
results.by.lev.master = list()
load(file=file.path(outputdir,"compare-methods-fourjump-123.Rdata"))
results.by.lev.master[1:3] = results.by.lev[1:3]
load(file=file.path(outputdir, "compare-methods-fourjump-45.Rdata"))
results.by.lev.master[4:5] = results.by.lev[1:2]
load(file=file.path(outputdir, "compare-methods-fourjump-67.Rdata"))
results.by.lev.master[6:7] = results.by.lev[1:2]
load(file=file.path(outputdir, "compare-methods-fourjump-89.Rdata"))
results.by.lev.master[8:9] = results.by.lev[1:2]
results.by.lev = results.by.lev.master

## Parse the results
myclean <- function(myresult){
    aa = lapply(myresult, function(a) do.call(rbind,a))
    return(aa)
}
mycleanresult = lapply(results.by.lev, myclean)

## Extract powers from it
levs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
names(mycleanresult) = levs
cond.pows = sapply(levs, function(mylev){
    pvs = mycleanresult[[toString(mylev)]][[1]][,1]
    pvs = pvs[!is.na(pvs)]
    mypow = sum(pvs<0.05/4)/length(pvs)
})

## Collect uncond pows
all.names = c("fl.rand", "fl.nonrand", "sbs.rand",
                 "sbs.nonrand", "wbs.rand", "wbs.nonrand",
                 "cbs.rand", "cbs.nonrand")[c(1,3,5,7)]

## cond.pows.by.method = sapply(all.names, function(methodname){
## uncond.pows.by.method = sapply(all.names, function(methodname){
recoveries.by.method = sapply(all.names, function(methodname){


    ## ## Investigating a few things:
    ## methodname="wbs.nonrand"
    ## methodname="fl.nonrand"
    ## par(mfrow=c(3,3))
    ## lapply(levs, function(mylev){
    ## pvs = mycleanresult[[toString(mylev)]][[methodname]][,"pvs"]
    ## locs = mycleanresult[[toString(mylev)]][[methodname]][,"locs"]
    ## hist(abs(locs)[!is.na(locs)])
    ## })

    methodname = "fl.rand"
    cond.pows = sapply(levs, function(mylev){
        print(mylev)
        pvs = mycleanresult[[toString(mylev)]][[methodname]][,"pvs"]
        pvs = pvs[!is.na(pvs)]
        mypow = sum(pvs<0.05/4)/length(pvs)
        print(mypow)
    })
    names(cond.pows) = levs
    return(cond.pows)


    ## uncond.pows = sapply(levs, function(mylev){
    ##     print(mylev)
    ##     pvs = mycleanresult[[toString(mylev)]][[methodname]][,"pvs"]
    ##     len = length(pvs)
    ##     pvs = pvs[!is.na(pvs)]
    ##     mypow = sum(pvs<0.05/4)/len
    ## })
    ## names(uncond.pows) = levs
    ## return(uncond.pows)

    ## recoveries = sapply(levs, function(mylev){
    ##     print(mylev)
    ##     len= nrow(mycleanresult[[toString(mylev)]][[methodname]])
    ##     mycleanresult[[toString(mylev)]][["sbs.nonrand"]]
    ##     pvs = mycleanresult[[toString(mylev)]][[methodname]][,"pvs"]
    ##     pvs = pvs[!is.na(pvs)]
    ##     ## mypow = sum(pvs<0.05/4)/length(pvs)
    ##     recovery = length(pvs)/(len)
    ## })
    ## names(recoveries) = levs
    ## return(recoveries)

})
cols = RColorBrewer::brewer.pal(4,"Set3")
lwd = rep(2,4)
lty = rep(1,4)
matplot(y=cond.pows.by.method, x=levs, type='l', col=cols, lwd=lwd, lty=lty)
matplot(y=recoveries.by.method, x=levs, type='l', col=cols, lwd=lwd, lty=lty)
legend("bottomright", col=cols, lwd=lwd, lty=lty, legend=all.names)



## Get detection ability
## lapply(results.by.lev[[2]]$wbs.rand, function(myresult){
##     any(is.na(unlist(myresult)))})
wbs.rand.detect.prop <- lapply(results.by.lev, function(my.result.by.lev){
    len = length(my.result.by.lev$wbs.rand)
    pvs = unlist(lapply(my.result.by.lev$wbs.rand, function(myresult){
        (myresult)[,"pvs"] }))
    pvs = pvs[!is.na(pvs)]
    recovery = length(pvs)/(4*len)
    return(recovery)
})

plot(pows, type='l')





## What does the onecompare() function produce NA's for?


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


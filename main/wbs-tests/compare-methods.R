## Synopsis: See randomized tests power compared to the nonrandomized version
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


jj = 3
whichlev.list = list(1:3, 4:6, 7:9)
whichlev = whichlev.list[[jj]]
levs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)[whichlev]
results.by.lev = list()
mc.cores = 8
nsims=c(seq(from=3000,to=1000,length=5), round(seq(from=600, to=300, length=4) ))[whichlev]
n=200 ## n=50
visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
print(levs)
for(ilev in 1:length(levs)){
    mylev = levs[ilev]
    nsim = nsims[ilev]
    print(mylev)
    results.by.lev[[ whichlev[ilev] ]] = onecompare(lev=mylev,
                                        nsim=nsim, meanfun=fourjump, visc=visc.fourjump,
                                        numSteps=4, bits=3000, mc.cores=mc.cores, n=n, numIS=100,
                                        max.numIS=3000)
    filename = paste0("compare-methods-fourjump-", paste0(whichlev.list[[jj]], collapse=""), ".Rdata")
    save(list=c("results.by.lev","levs","nsim", "n"), file=file.path(outputdir, filename))
    print(filename)
}


## Why is fl.nonrand failing? ## TODOcheck and erase
myname = "fl.nonrand"
a = mclapply(1:nsim, function(isim) {
    printprogress(isim,nsim);
    dosim_compare(type=myname, n=n, lev=lev, numIS=numIS,
                  meanfun=meanfun, visc=visc, numSteps=numSteps, bits=bits)
}, mc.cores=mc.cores)



## Aggregate the results
outputdir = "../output"
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
## load(file=file.path(outputdir,"compare-methods-fourjump-123.Rdata")) ## Temporary, for flplus 3
mycleanresult = lapply(results.by.lev, myclean)
levs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
names(mycleanresult) = levs


## Collect uncond pows
all.names = c("fl.rand", "fl.rand.plus",  "fl.nonrand", "sbs.rand",
                 "sbs.nonrand", "wbs.rand", "wbs.nonrand",
                 "cbs.rand", "cbs.nonrand")#[c(2,4)]#[c(1,3,5,7)]

cond.pows.by.method = sapply(all.names, function(methodname){
    cond.pows = sapply(levs, function(mylev){
        pvs = mycleanresult[[toString(mylev)]][[methodname]][,"pvs"]
        pvs = pvs[!is.na(pvs)]
        mypow = sum(pvs<0.05/4)/length(pvs)
    })
    names(cond.pows) = levs
    return(cond.pows)
})

uncond.pows.by.method = sapply(all.names, function(methodname){
    uncond.pows = sapply(levs, function(mylev){
        pvs = mycleanresult[[toString(mylev)]][[methodname]][,"pvs"]
        len = nsims[toString(mylev)]
        ## len=300
        pvs = pvs[!is.na(pvs)]
        mypow = sum(pvs<0.05/4)/(len*4)
    })
    print(uncond.pows)
    names(uncond.pows) = levs
    return(uncond.pows)
})

recoveries.by.method = sapply(all.names, function(methodname){
    recoveries = sapply(levs, function(mylev){
        print(mylev)
        len = nsims[toString(mylev)]
        mycleanresult[[toString(mylev)]][["sbs.nonrand"]]
        pvs = mycleanresult[[toString(mylev)]][[methodname]][,"pvs"]
        pvs = pvs[!is.na(pvs)]
        recovery = length(pvs)/(4*len)
    })
    names(recoveries) = levs
    return(recoveries)
})

par(mfrow=c(3,1))
cols = RColorBrewer::brewer.pal(4,"Set2")
lwd = rep(2,4)
lty = rep(1,4)
ylim = c(0,1)
xlab = expression(delta)
ylab = ""

pdf(file=file.path(outputdir,"compare-methods-cond.pdf"), width=5, height=5)
matplot(y=cond.pows.by.method, x=levs, type='l', col=cols, lwd=lwd, lty=lty, main = "Conditional Power", axes=FALSE, ylim=ylim, ylab=ylab, xlab=xlab)
axis(2);axis(1)
legend("bottomright", col=cols, lwd=lwd, lty=lty, legend=all.names)
graphics.off()


pdf(file=file.path(outputdir,"compare-methods-uncond.pdf"), width=5, height=5)
matplot(y=uncond.pows.by.method, x=levs, type='l', col=cols, lwd=lwd, lty=lty, main = "Unconditional Power",axes=FALSE, ylim=ylim, ylab=ylab, xlab=xlab)
axis(2);axis(1)
legend("bottomright", col=cols, lwd=lwd, lty=lty, legend=all.names)
graphics.off()

pdf(file=file.path(outputdir,"compare-methods-recov.pdf"), width=5, height=5)
matplot(y=recoveries.by.method, x=levs, type='l', col=cols, lwd=lwd, lty=lty, main = "Detection Probability", axes=FALSE, ylim=ylim, ylab=ylab, xlab=xlab)
axis(2);axis(1)
legend("bottomright", col=cols, lwd=lwd, lty=lty, legend=all.names)
graphics.off()


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


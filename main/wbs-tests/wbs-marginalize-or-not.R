## Synopsis: Want to see if marginalization of WBS inference is worth it, in
## terms of (conditional) power! Taking results from compare-methods

outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
cols = RColorBrewer::brewer.pal(2, "Set1")
w=10; h=7
## loc = c(40,80,120,160)
## pdf(file=file.path("~/Desktop", paste0("wbs-around-any-loc.pdf" )), width=w,height=h)
loc = 160
## pdf(file=file.path("~/Desktop", paste0("wbs-around-", loc, ".pdf" )), width=w,height=h)
par(mfrow=c(2,4))
levs.master = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
for(jj in 1:4){
    whichlev.list = list(1:3, 4:5,6:7, 8:9)
    whichlev.list2 = list(2:3, 4:5,6:7, 8:9)
    filename = paste0("compare-methods-fourjump-", paste0(whichlev.list[[jj]], collapse=""), ".Rdata")
    load(file=file.path(outputdir, filename))
    for(kk in whichlev.list2[[jj]]){
        print(kk)
        rand.results = (results.by.lev[[kk]])$wbs.rand
        nonrand.results = (results.by.lev[[kk]])$wbs.nonrand

        loc.vicinity = as.numeric(sapply(loc, function(myloc) myloc + c(-1,0,1)))

        pvs.rand = unlist(sapply(rand.results, function(a) a$pvs))
        locs.rand = unlist(sapply(rand.results, function(a) a$locs))
        pvs.rand = pvs.rand[abs(locs.rand) %in% loc.vicinity]
        pvs.rand[is.nan(pvs.rand)] = 0
        pvs.rand = pvs.rand[which(!is.na(pvs.rand))]

        pvs.nonrand = unlist(sapply(nonrand.results, function(a) a$pvs))
        locs.nonrand = unlist(sapply(nonrand.results, function(a) a$locs))
        pvs.nonrand = pvs.nonrand[abs(locs.nonrand) %in% loc.vicinity]
        pvs.nonrand = pvs.nonrand[which(!is.na(pvs.nonrand))]

        qqunif(list(marginalized=pvs.rand, not.marginalized=pvs.nonrand),cols=cols)
        lev = levs.master[kk]
        title(main=bquote(delta==.(lev)))
    }
}
graphics.off()


## Also plot the power

## Aggregate the results from the computers
outputdir = "../output"
results.by.lev.master=list()
for(jj in 1:4){
    filename = paste0("compare-methods-fourjump-", paste0(whichlev.list[[jj]], collapse=""), ".Rdata")
    load(file=file.path(outputdir, filename))
    results.by.lev.master[whichlev.list[[jj]] ] = results.by.lev[whichlev.list[[jj]]]
}
nsims=c(seq(from=3000,to=1000,length=5), round(seq(from=600, to=300, length=4) ))
all.names.list = list( c( "wbs.rand", "wbs.nonrand"),  c( "sbs.rand", "sbs.nonrand"))

## Make plot
w=10
h=5
pdf(file=file.path(outputdir, "cond-power-comparison.pdf"), width=w, height=h)
par(mfrow=c(1,2))
cond.pows.by.method.list = list()
cols.list = RColorBrewer::brewer.pal(2, "Set2")
levs.master = levs.master[-1]
for(iplot in 1:2){
    all.names = all.names.list[[iplot]]

    ## Parse the results
    myclean <- function(myresult){  return(lapply(myresult, function(a) do.call(rbind,a))) }
    mycleanresult = lapply(results.by.lev.master, myclean)
    names(mycleanresult) = levs.master


    ## Collect uncond pows
    cond.pows.by.method = sapply(all.names, function(methodname){
        cond.pows = sapply(levs.master, function(mylev){
            pvs = mycleanresult[[toString(mylev)]][[methodname]][,"pvs"]
            pvs = pvs[!is.na(pvs)]
            mypow = sum(pvs<0.05/4)/length(pvs)
        })
        names(cond.pows) = levs.master
        return(cond.pows)
    })
    cond.pows.by.method.list[[iplot]] = cond.pows.by.method

    cols = cols.list[[iplot]]
    matplot(y=cond.pows.by.method,x=levs.master, lwd=2, type='o', lty=1, ylim=c(0,1), pch=c(1,2), col=cols, cex=1)

    if(iplot==1){
        legend("bottomright", lwd=c(3,3),
               legend=c("Marginalized WBS", "Nonmarginalized WBS"), lty=c(1,1),
               pch = c(1,2), col = cols, cex=1)
        title(main="Marginalized vs Nonmarginalized WBS")
    } else {
        legend("bottomright", lwd=c(3,3),
               legend=c("Noise-added binseg", "Plain binseg"), lty=c(1,1),
               pch = c(1,2), col = cols, cex=1)
        title(main="Noise-added BS vs BS")
    }
}
graphics.off()

w=h=5
pdf(file=file.path(outputdir, "cond-power-comparison-combined.pdf"), width=w, height=h)
par(mfrow=c(1,1))
par(mar=2,2,2,2)
cond.pows = do.call(cbind, cond.pows.by.method.list)
matplot(y=cond.pows, x=levs.master, lwd=2, type='o', lty=1, ylim=c(0,1), pch=c(1,2,1,2), col=cols.list[c(1,1,2,2)], cex=1)
legend("bottomright", lwd=2,
       legend=c("Marg. WBS", "Nonmarg. WBS", "Noise-added binseg", "Plain binseg"), lty=c(1,1,1,1),
       pch = c(1,2,1,2), col = rep(cols.list, each=2), cex=1)
graphics.off()

uncond.pows.by.method = sapply(all.names, function(methodname){
    uncond.pows = sapply(levs.master, function(mylev){
        pvs = mycleanresult[[toString(mylev)]][[methodname]][,"pvs"]
        len = nsims[toString(mylev)]
        ## len=300
        pvs = pvs[!is.na(pvs)]
        mypow = sum(pvs<0.05/4)/(len*4)
    })
    print(uncond.pows)
    names(uncond.pows) = levs.master
    return(uncond.pows)
})

recoveries.by.method = sapply(all.names, function(methodname){
    recoveries = sapply(levs.master, function(mylev){
        print(mylev)
        len = nsims[toString(mylev)]
        mycleanresult[[toString(mylev)]][["sbs.nonrand"]]
        pvs = mycleanresult[[toString(mylev)]][[methodname]][,"pvs"]
        pvs = pvs[!is.na(pvs)]
        recovery = length(pvs)/(4*len)
    })
    names(recoveries) = levs.master
    return(recoveries)
})


## Plot the results
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


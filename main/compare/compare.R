## Synopsis: See each method's power 
library(genlassoinf)
outputdir = "../output"
source("../main/compare/sim-helpers.R")
onecompare <- function(lev=0, nsim=1000, mc.cores=8, meanfun=onejump, visc=NULL,
                       numSteps=1, bits=50, n=60, numIS=200,
                       max.numIS=1000){

    all.names = type=c("fl","fl.noisy","fl.noisy.plus",
                       "sbs.noisy", "sbs",
                       "wbs.marg", "wbs",
                       "cbs", "cbs.noisy")

    ## all.results = lapply(all.names, function(myname){
    ##     print(myname)
    all.results = list()
    for(iname in 1:length(all.names)){
        myname = all.names[iname]
        myresult = NULL
        tryCatch({
            myresult =  mclapply(1:nsim, function(isim) {
                printprogress(isim,nsim);
                dosim_compare(type=myname, n=n, lev=lev, numIS=numIS,
                              meanfun=meanfun, visc=visc, numSteps=numSteps,
                              bits=bits, max.numIS=max.numIS)
            }, mc.cores=mc.cores)
            
        }, error = function(err) {
            err$message = paste(err$message,"\n(", name, " resulted in an error.)",sep="")
            warning(err)
        })
        cat(fill=TRUE) 
        all.results[[myname]] = myresult
    }
    names(all.results) = all.names
    return(all.results)
}

myresult = onecompare(lev=mylev, nsim=nsim, meanfun=fourjump, visc=visc.fourjump,
                      numSteps=4, bits=3000, mc.cores=mc.cores, n=n, numIS=100,
                      max.numIS=3000)


## Run the actual simulations
levs = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4)[whichlev]
nsims = c(3000,3000,3000, seq(from=3000,to=1000,length=5),
            round(seq(from=600, to=300, length=4) ))[whichlev]
n = 200
visc.fourjump = unlist(lapply(c(1,2,3,4) * (n/5), function(cp){ cp+c(-1,0,1)}))
results.by.lev = list()
for(ilev in 1:length(levs)){
    cat(fill=TRUE)
    mylev = levs[ilev];  nsim = nsims[ilev]
    printprogress(mylev, levs, "levels (jump size)")
    thisresult = onecompare(lev=mylev, nsim=nsim, meanfun=fourjump, visc=visc.fourjump,
                            numSteps=4, bits=3000, mc.cores=mc.cores, n=n, numIS=100,
                            max.numIS=3000)
    results.by.lev[[ whichlev[ilev] ]] = thisresult
    filename = paste0("compare-fourjump-", ilev, ".Rdata")
    save(list=c("results.by.lev","levs","nsim", "n"), file=file.path(outputdir, filename))
    print(filename)
}

## Trying out one method at a time, to make sure they work. 
n = 20
visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
mylev=1
mc.cores=4
lev=mylev; nsim=nsim; meanfun=fourjump; visc=visc.fourjump;
numSteps=4; bits=3000; mc.cores=mc.cores; n=n; numIS=100;
max.numIS=3000
myname = "fl.noisy"
nsim=100
for(isim in 1:nsim){
    printprogress(isim, nsim)
    cat(fill=TRUE)
    a = dosim_compare(type=myname, n=n, lev=lev, numIS=numIS, meanfun=meanfun,
                      visc=visc, numSteps=numSteps, bits=bits, max.numIS=max.numIS,
                      verbose=FALSE)
    cat(fill=TRUE)
}

mc.cores=4
nsim=4
n=20
thisresult = onecompare(lev=mylev, nsim=nsim, meanfun=fourjump,
                        visc=visc.fourjump, numSteps=4, bits=3000,
                        mc.cores=mc.cores, n=n, numIS=100, max.numIS=3000)



## What do I need to check? that the IS is all done to some sampling precision.



## Aggregate the results from the computers
outputdir = "../output"
for(jj in 1:3){
    jj=1
    filename = paste0("compare-methods-fourjump-", paste0(whichlev.list[[jj]], collapse=""), ".Rdata")
    load(file=file.path(outputdir, filename))
    results.by.lev.master[jj] = results.by.lev[jj]
}


## Parse the results
myclean <- function(myresult){  return(lapply(myresult, function(a) do.call(rbind,a))) }
mycleanresult = lapply(results.by.lev, myclean)
names(mycleanresult) = levs


## Collect uncond pows
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




#####################
### Other stuff #####
#####################

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

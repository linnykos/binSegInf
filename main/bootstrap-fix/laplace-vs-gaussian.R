## Synopsis: Checking out systematically where bootstrap go wrong. Change the file name later!
library(smoothmest)
la("~/repos/genlassoinf/genlassoinf/")
source("../main/bootstrap-fix/helper.R")
outputdir = "../output"

## Run two types of inferences
nsim = 5000
n = 200
lev = 1
mc.cores = 6
max.numSteps = 10
pvs.gaus.by.step =list()
pvs.laplace.by.step =list()
start.time = Sys.time()
for(ii in 1:max.numSteps){
    printprogress(ii,max.numSteps,"numSteps",fill=TRUE, start.time=start.time)
    pvs.gaus = unlist(mclapply(1:nsim, function(isim){
        dosim(n=n,numSteps=ii,type="gaus", lev=lev)}, mc.cores=mc.cores))
    pvs.laplace = unlist(mclapply(1:nsim, function(isim){
        dosim(n=n,numSteps=ii,type="laplace", lev=lev)}, mc.cores=mc.cores))
    pvs.gaus.by.step[[ii]] = pvs.gaus
    pvs.laplace.by.step[[ii]] = pvs.laplace
}
filename = "bootstrap-laplace-by-step-lev0.Rdata"
filename = "bootstrap-laplace-by-step-lev1.Rdata"
save(pvs.laplace.by.step, pvs.gaus.by.step, file=file.path(outputdir, filename))


## Plot the results in QQ plots
load(file=file.path(outputdir, filename))
plotfilename = makeplotfilename(filename)
## pdf(file=file.path(outputdir, plotfilename, width=15, height=10)
jpeg(file=file.path(outputdir, plotfilename), width=2000, height=1000)
par(mfrow=c(3, 4))
for(ii in 1:max.numSteps){
    qqunif(list(lapl=pvs.laplace.by.step[[ii]],
                gaus=pvs.gaus.by.step[[ii]]), col=c("black","red"))
    title(main=paste0("step", ii))
}
qqunif(pvs.laplace.by.step, cols=1:max.numSteps); legend("bottomright", col=1:5, pch=16, legend=1:5)
title(main="All steps")
type1errors = sapply(pvs.laplace.by.step, function(pvs)sum(pvs<0.05)/length(pvs))
plot(type1errors, type='l', ylim = c(0,0.2), lwd=2)
title(main="Proportion of p-values under 0.05")
abline(h=0.05, lwd=2, lty=2)
graphics.off()


## Examine IC stopping inference
nsim = 5000
n = 200
mc.cores = 6
max.numSteps = 10
lev = 0
start.time = Sys.time()
pvs.laplace.test = unlist(mclapply(1:nsim, function(isim){
    printprogress(isim,nsim, start.time=start.time)
    dosim_icstop(n=n,numSteps=10, type="laplace", lev=lev)}, mc.cores=mc.cores))
pvs.gaus.test = unlist(mclapply(1:nsim, function(isim){
    printprogress(isim,nsim, start.time=start.time)
    dosim_icstop(n=n,numSteps=10, type="gaus", lev=lev)}, mc.cores=mc.cores))
filename = "bootstrap-laplace-icstop-lev0.Rdata"
## filename = "bootstrap-laplace-icstop-lev1.Rdata"
save(pvs.gaus.test, pvs.laplace.test, file=file.path(outputdir, filename))


## Load and plot
load(file.path(outputdir, filename))
plotfilename = makeplotfilename(filename)
jpeg(file=file.path(outputdir, plotfilename), width=500, height=500)
qqunif(list(lapl=pvs.laplace.test,
            gaus=pvs.gaus.test), col=c("black","red"))
graphics.off()


## Example of data
jpeg(file=file.path(outputdir, "onejump-data-example.jpg"), width=500, height=500)
lev = 1
n = 200
mn = c(rep(0, n/2), rep(lev, n/2))
y = mn + rnorm(n,0,1)
plot(y, pch=16, col='grey50')
lines(mn, lwd=2, col='red')
graphics.off()

## Why does fused lasso do so well? Why did binseg have so many null hits after
## IC before?
mc.cores=4
nsim=1000
lev=1
n=200
start.time = Sys.time()
stoptimes.gaus.binseg = unlist(mclapply(1:nsim, function(isim){
    printprogress(isim,nsim, start.time=start.time)
    dosim_icstop(n=n,numSteps=10, type="gaus", lev=lev, model="binseg")}, mc.cores=mc.cores))

stoptimes.gaus.fusedlasso = unlist(mclapply(1:nsim, function(isim){
    printprogress(isim,nsim, start.time=start.time)
    dosim_icstop(n=n,numSteps=10, type="gaus", lev=lev, model="fusedlasso")}, mc.cores=mc.cores))

## This significantly different?
## This is interesting; why would this be happening?
plot(table(stoptimes.gaus.fusedlasso))
lines(y=table(stoptimes.gaus.binseg),x=(1:8),pch=16)


## After checking this, do the fused lasso 
nsim = 5000
n = 200
mc.cores = 7
lev = 0
start.time = Sys.time()
pvs.laplace.test = unlist(mclapply(1:nsim, function(isim){
    printprogress(isim,nsim, start.time=start.time)
    dosim_icstop(n=n,numSteps=10, type="laplace", lev=lev, model="fusedlasso")}, mc.cores=mc.cores))
pvs.gaus.test = unlist(mclapply(1:nsim, function(isim){
    printprogress(isim,nsim, start.time=start.time)
    dosim_icstop(n=n,numSteps=10, type="gaus", lev=lev, model="fusedlasso")}, mc.cores=mc.cores))
## filename = "bootstrap-laplace-icstop-fusedlasso-lev0.Rdata"
filename = "bootstrap-laplace-icstop-fusedlasso-lev1.Rdata"
save(pvs.gaus.test, pvs.laplace.test, file=file.path(outputdir, filename))


## Load and plot
load(file=file.path(outputdir, filename))
plotfilename = makeplotfilename(filename)
## pdf(file=file.path(outputdir, plotfilename, width=15, height=10)
jpeg(file=file.path(outputdir, plotfilename), width=500, height=500)
qqunif(list(lapl=pvs.laplace.test,
            gaus=pvs.gaus.test), col=c("black","red"))
graphics.off()



## On fused lasso, for 1-10 steps, run two types of inferences:
nsim = 5000
n = 200
lev = 0
mc.cores = 6
max.numSteps = 10
pvs.gaus.by.step =list()
pvs.laplace.by.step =list()
start.time = Sys.time()
for(ii in 1:max.numSteps){
    printprogress(ii,max.numSteps,"numSteps",fill=TRUE, start.time=start.time)
    pvs.gaus = unlist(mclapply(1:nsim, function(isim){
        dosim(n=n,numSteps=ii,type="gaus", lev=lev, model="fusedlasso")}, mc.cores=mc.cores))
    pvs.laplace = unlist(mclapply(1:nsim, function(isim){
        dosim(n=n,numSteps=ii,type="laplace", lev=lev, model="fusedlasso")}, mc.cores=mc.cores))
    pvs.gaus.by.step[[ii]] = pvs.gaus
    pvs.laplace.by.step[[ii]] = pvs.laplace
}
filename = "bootstrap-laplace-by-step-fusedlasso-lev0.Rdata"
filename = "bootstrap-laplace-by-step-fusedlasso-lev1.Rdata"
save(pvs.laplace.by.step, pvs.gaus.by.step, file=file.path(outputdir, filename))


## Plot the results in QQ plots
load(file=file.path(outputdir, filename))
plotfilename = makeplotfilename(filename)
## pdf(file=file.path(outputdir, plotfilename, width=15, height=10)
jpeg(file=file.path(outputdir, plotfilename), width=2000, height=1000)
par(mfrow=c(3, 4))
for(ii in 1:max.numSteps){
    qqunif(list(lapl=pvs.laplace.by.step[[ii]],
                gaus=pvs.gaus.by.step[[ii]]), col=c("black","red"))
    title(main=paste0("step", ii))
}
qqunif(pvs.laplace.by.step, cols=1:max.numSteps); legend("bottomright", col=1:5, pch=16, legend=1:5)
title(main="All steps")
type1errors = sapply(pvs.laplace.by.step, function(pvs)sum(pvs<0.05)/length(pvs))
plot(type1errors, type='l', ylim = c(0,0.2), lwd=2)
title(main="Proportion of p-values under 0.05")
abline(h=0.05, lwd=2, lty=2)
graphics.off()


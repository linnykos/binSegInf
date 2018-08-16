## Synopsis: Checking out systematically where bootstrap go wrong. Change the file name later!
library(smoothmest)
outputdir = "../output"

## Helper:
dosim <- function(n=200, numSteps=1, lev=0, type=c("gaus", "laplace", "bootstrap"), retain.null=FALSE,
                  icstop=TRUE){
    mn = c(rep(0, n/2), rep(lev, n/2))
    if(type=="laplace") y = mn + rdoublex(n,lambda=1/sqrt(2))
    if(type=="gaus") y = mn + rnorm(n,0,1)
    h = binSeg_fixedSteps(y, numSteps=numSteps)
    h.poly = polyhedra(h, numSteps=numSteps)
    vlist = make_all_segment_contrasts_from_cp(cp=h$cp, cp.sign=h$cp.sign, n=n)
    if(retain.null){
        retain = which(sapply(vlist, function(v){all.equal(sum(v*mn),0)==TRUE}))
        if(length(retain)==0) return(NULL)
    } else {
        retain = 1:length(vlist)
    }
    vlist = vlist[retain]
    pvs = sapply(vlist, function(v){poly.pval2(y=y, v=v, poly=h.poly, sigma=1)$pv})
    return(pvs)
}

## Run two types of inferences
nsim = 5000
n = 200
lev = 1
mc.cores = 6
max.numSteps = 10
retain.null = TRUE
pvs.gaus.by.step =list()
pvs.laplace.by.step =list()
start.time = Sys.time()
for(ii in 1:max.numSteps){
    printprogress(ii,max.numSteps,"numSteps",fill=TRUE, start.time=start.time)
    pvs.gaus = unlist(mclapply(1:nsim, function(isim){
        dosim(n=n,numSteps=ii,type="gaus", lev=lev, retain.null=retain.null)}, mc.cores=mc.cores))
    pvs.laplace = unlist(mclapply(1:nsim, function(isim){
        dosim(n=n,numSteps=ii,type="laplace", lev=lev, retain.null=retain.null)}, mc.cores=mc.cores))
    pvs.gaus.by.step[[ii]] = pvs.gaus
    pvs.laplace.by.step[[ii]] = pvs.laplace
}
## filename = "bootstrap-laplace-by-step-lev0.Rdata"
filename = "bootstrap-laplace-by-step-lev1.Rdata"
save(pvs.laplace.by.step, pvs.gaus.by.step, file=file.path(outputdir, filename))


## Plot the results in QQ plots
load(file=file.path(outputdir, filename))
pdf(file=file.path("~/Desktop", "laplace.pdf"), width=10, height=8)
par(mfrow=c(2,4))
for(ii in 1:max.numSteps){
    qqunif(list(lapl=pvs.laplace.by.step[[ii]],
                gaus=pvs.gaus.by.step[[ii]]), col=c("black","red"))
}
qqunif(pvs.laplace.by.step, cols=1:max.numSteps); legend("bottomright", col=1:5, pch=16, legend=1:5)
type1errors = sapply(pvs.laplace.by.step, function(pvs)sum(pvs<0.05)/length(pvs))
plot(type1errors, type='l', ylim = c(0,0.5), lwd=2)
abline(h=0.03, lwd=2, lty=2)
graphics.off()




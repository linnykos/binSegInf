## Synopsis: Try to make a better, improved segment test by using the winning
## intervals as segments

## Sample settings
set.seed(1)
sigma=1
lev=3
n=10
meanfun = onejump
numSteps=1
numIntervals=n
randomized=FALSE
visc = n/2 + c((-1),0,+1)
locs = visc

n=60
nsims=c(3000,700,500,250)
## numSteps=3
mc.cores=1
visc = (n/2+((-1):1))
## results = lapply(levs, dosim, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump, mc.cores=mc.cores)
source("../main/wbs-tests/sim-helpers.R")
results = Map(function(lev,nsim)dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps,
                                      randomized=randomized,numIS=100,meanfun=onejump,
                                      mc.cores=mc.cores, locs=visc, better.segment=TRUE), levs, nsims)



## On fourjump example
source("../main/wbs-tests/sim-helpers.R")
whichlevs = 3:5
levs = c(0,1,2,3,4)[whichlevs]
n=200
mc.cores=8
meanfun = fourjump
numSteps = 4
randomized=TRUE
## visc = n/5 * (1:4)
visc = n/5 + c((-2):2)
mc.cores=8
nsims = seq.int(from=1000, to=200, length=5)[whichlevs]
for(ilev in 1:length(whichlevs)){
    lev = levs[ilev]
    nsim = nsims[ilev]
    print(lev)
    results = dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps, randomized=randomized,
                numIS=100, meanfun=meanfun, mc.cores=mc.cores, locs=visc,
                better.segment=TRUE, improve.nomass.problem = TRUE,
                min.num.things=20)
    results.orig = dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps, randomized=randomized,numIS=100,
                         meanfun=meanfun, mc.cores=mc.cores, locs=visc, better.segment=FALSE,
                         improve.nomass.problem = TRUE, min.num.things=20)
    filename = paste0("better-segment-fourjump-lev", lev, ".Rdata")
    save(list=c("results", "results.orig"), file=file.path(outputdir, filename))
}

## Load and plot better vs original segment test.
load(file=file.path(outputdir, filename))
qqunif(results$pv)
a = qqunif(results.orig$pv, plot.it=FALSE)
points(a,col='red')
legend("bottomright", legend=c("orig", "improved"), col = c("red", "black"), pch=c(16,16))





## ## On onejump example
## source("../main/wbs-tests/sim-helpers.R")
## outputdir = "../output"
## lev=1
## n=10
## mc.cores=4
## nsim=300
## meanfun = onejump
## numSteps=1
## randomized=TRUE
## ## visc = n/5 * (1:4)
## visc = n/2 + (c((-1):1))
## numIS=100
## improve.nomass.problem = TRUE
## results = dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps, randomized=randomized,numIS=numIS,
##                 meanfun=meanfun, mc.cores=mc.cores, locs=visc, better.segment=TRUE,
##                 improve.nomass.problem= improve.nomass.problem )
## results.orig = dosim(lev=lev,n=n,nsim=nsim,numSteps=numSteps, randomized=randomized,numIS=numIS,
##                      meanfun=meanfun, mc.cores=mc.cores, locs=visc, better.segment=FALSE,
##                      improve.nomass.problem =improve.nomass.problem )


## filename = "better-segment-onejump.Rdata"
## save(list=c("results.orig"), file=file.path(outputdir, filename))
## load(file=file.path(outputdir, filename))
## qqunif(results$pv)
## a = qqunif(results.orig$pv, plot.it=FALSE)
## points(a,col='blue')
## legend("bottomright", legend=c("orig", "improved"), col = c("blue", "black"), pch=c(16,16))


## Example figure
set.seed(0)
sigma = 1
n=200

## Trying a few examples out
set.seed(3)
y0 = fourjump(lev=1,n=n) + rnorm(n,0,sigma)
g = wildBinSeg_fixedSteps(y=y0, numSteps=3, numIntervals=n, inference.type="none")
plot(y0, ylim = ylim,axes=F, xlim=xlim, xlab = xlab, ylab = "", pch=pch, col=pcol);
abline(v=g$cp)
abline(v=g$results[which(g$results[,"max.b"]==85),c("max.s", "max.e")], col='blue')

nsim=100
jumps=matrix(nrow=nsim, ncol=4)
closest.ends = ends=matrix(nrow=nsim,ncol=2)
for(isim in 1:nsim){
    y0 = fourjump(lev=1,n=n) + rnorm(n,0,sigma)
    g = wildBinSeg_fixedSteps(y=y0, numSteps=4, numIntervals=n, inference.type="none")
    jumps[isim,] = g$cp
    if(any(g$cp %in% c(79,80,81))){
        ii = which(g$cp %in% c(79,80,81))
        ends[isim,] = g$results[ii,c("max.s", "max.e")]
        which.closest = order(abs(g$cp - g$cp[ii]))[2:3]
        closest.ends[isim,] = sort(g$cp[which.closest])
    }
}

cbind(ends,closest.ends)
par(mfrow=c(2,2))
hist(ends[,1], xlim=c(0,n))
abline(v=c(40,120), col='red')
hist(ends[,2], xlim=c(0,n))
abline(v=c(40,120), col='red')
hist(closest.ends[,1], xlim=c(0,n))
abline(v=c(40,120), col='red')
hist(closest.ends[,2], xlim=c(0,n))
abline(v=c(40,120), col='red')






beta0 = fourjump(lev=1,n=n)
brk = 40
brk.short = 20

v.orig = c(.3*rep(-0.7,brk) , .3*rep(+0.7,2*brk) )
x.contrasts.orig = c(1:120)
v.improved = c(rep(NA,40-brk.short), c(.3*rep(-0.7,brk.short) , .3*rep(+0.7,brk.short) )*2, rep(NA,40-brk.short))
x.contrasts.improved = c(1:80)

xlab = "Location"
w = 5; h = 5
pch = 16; lwd = 2
pcol = "gray50"
ylim = c(-5,5)
mar = c(4.5,4.5,0.5,0.5)
xlim = c(0,n+10)

xticks = c(0,2,4,6)*10
let = c("A","B")
ltys.sig = c(2,2,1)
lwd.sig = 2
pch.dat = 16
pcol.dat = "grey50"
pch.contrast = 17
lty.contrast = 2
lcol.sig = 'red'
pcol.spike=3
pcol.segment=4
pcols.delta =   pcols.oneoff = RColorBrewer::brewer.pal(n=3,name="Set2")
pch.spike = 15
pch.segment = 17
cex.contrast = 1.2

Plot single example of data and contrast (top leftmost figure)
pdf(file.path(outputdir,"better-contrast-data.pdf"), width=5,height=5)

## Run Wild binary segmentation for four steps
g = wildBinSeg_fixedSteps(y=y0, numSteps=4, numIntervals=n, inference.type="none")

par(mar=c(4.1,3.1,3.1,1.1))
 plot(y0, ylim = ylim,axes=F, xlim=xlim, xlab = xlab, ylab = "", pch=pch, col=pcol);
 axis(1, at = xticks, labels = xticks); axis(2)
 ## for(ii in 1:3) lines(beta0[[ii]],col="red",lty=ltys.sig[ii], lwd=lwd.sig)
 ## for(ii in 0:2) text(x=65,y=ii, label = bquote(delta==.(ii)))
lines(beta0,col="red",lty=ltys.sig[3], lwd=lwd.sig)
 points(v.orig~x.contrasts.orig, pch = pch.spike, col = 3)
 points(v.improved~x.contrasts.improved, pch = pch.segment, col = 4)

 abline(h = mean(v.segment,na.rm=T), col = 'lightgrey')
 legend("topleft", pch=c(pch.dat,NA,pch.spike,pch.segment),
        lty=c(NA,1,NA,NA), lwd=c(NA,2,NA,NA),
        col = c(pcol.dat, lcol.sig, pcol.spike,pcol.segment),
        pt.cex = c(cex.contrast, NA, cex.contrast, cex.contrast),
        legend=c("Data", "Mean","Spike contrast", "Segment contrast"))
title(main=expression("Data example"))
## graphics.off()


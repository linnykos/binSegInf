source('../main/sat-vs-sel/sim-helper.R')

## Generate a toy dataset
n = 12
twojump <- function(lev,n){c(rep(0,n/3),rep(lev,n/3), rep(0,n/3))}
lev=1
sigma=1


## Plot settings
mar = c(4.5,4.5,2.5,0.5)
w=h=5
set.seed(0)
pch.dat=16
pcol.dat="grey50"
lcol.sig="red"
lwd.sig = 2
lev1 = 0
lev2 = 3
n = 20
sigma = 1
ylab = ""
xlab = "Coordinate"
lty.sig = 1
xlim = c(0,24)
ylim = c(-2,4)

outputdir = "../output"
pdf(file=file.path(outputdir,"sat-vs-sel-data.pdf"),width=w,height=h)
par(mar=mar)


## Generate data + Fit WBS
set.seed(1)
y = mn + rnorm(n, 0, sigma)
numSteps=2
## g = binSeg_fixedSteps(y, numSteps=numSteps)
mn = twojump(lev,n)
g = wildBinSeg_fixedSteps(y, numSteps=numSteps, numIntervals=n)
poly = polyhedra(obj=g$gamma, u=g$u)

plot(y, col = pcol.dat, axes=F,ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,pch=pch.dat)
axis(1);axis(2)
lines(mn, col=lcol.sig,lwd= lwd.sig)
legend("topleft",lty=c(lty.sig,NA), lwd=c(lwd.sig,NA),
        col = c(lcol.sig,pcol.dat), pch = c(NA,pch.dat),
        legend = c("Signal","Data"))
title(main="Example of data and signal \n(n=12)", line=1)
abline(v=10,col='lightgrey',lty=3)
graphics.off()

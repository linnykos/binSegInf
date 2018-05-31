## Synopsis: Renewed comparison of Vlo of one-step fused lasso and binary
## segmentation, focusing on different regions
source("helpers.R")
outputdir = "../figures"
la("~/repos/binSegInf/binSegInf")
la("~/repos/genlassoinf/genlassoinf/")

## Generate date
lev=100
mn = c(rep(0,5), rep(lev, 30), rep(0,5), rep(lev, 30))
set.seed(3)
y = mn + rnorm(length(mn), 0, 1)
D = makeDmat(length(y),order=0)
out = dualpathSvd4(y, D=D, maxsteps=4)


## Run simulation
lapl <- function(n){ rexp(n,rate=sqrt(2)) * sample(c(-1,1),n,replace=TRUE)}
nrep = 20000##500000
out.lapl = sim(n=100, nrep=nrep, err.fun=lapl, seed=NULL)
out.gaus = sim(n=100, nrep=nrep, err.fun=rnorm, seed=NULL)
## ## outputdir = "."
## save(out.lapl, out.gaus, file=file.path(outputdir, "vlo.Rdata")
## load(file=file.path(".", "vlo.Rdata"))


mean.vlo.gaus.bs = sapply(1:99, function(loc){
    mean(out.gaus$vlo.bs[which(out.gaus$j.bs==loc)])
})
mean.vlo.gaus.fl = sapply(1:99, function(loc){
    mean(out.gaus$vlo.fl[which(out.gaus$j.fl==loc)])
})
mean.vlo.lapl.bs = sapply(1:99, function(loc){
    mean(out.lapl$vlo.bs[which(out.lapl$j.bs==loc)])
})
mean.vlo.lapl.fl = sapply(1:99, function(loc){
    mean(out.lapl$vlo.fl[which(out.lapl$j.fl==loc)])
})

## makejpg(outputdir, "vlo-gaus.jpg")
## main = "Vlo, Gaussian"
## ylim = c(0,5)
## plot(mean.vlo.gaus.bs, type='o', pch=16, ylim=ylim, ylab = "location",
##      xlab = "", main=main, lwd=2)
## lines(mean.vlo.gaus.fl, type='o', pch=16, ylim=ylim, ylab = "location",
##      xlab = "", main=main, col='red', lwd=2)
## graphics.off()

makejpg(outputdir, "vlo-by-loc.jpg", 700,900)
main = "Vlo by location"
ylim = c(0,5)
cols=c(1,2)
plot(mean.vlo.lapl.bs, type='l', pch=16, ylim=ylim, ylab = "location",
     xlab = "", main=main, lwd=1)
lines(mean.vlo.lapl.fl, type='l', pch=16, ylim=ylim, ylab = "location",
     xlab = "", main=main, col='red', lwd=1)
lines(mean.vlo.gaus.bs, type='l', pch=16, ylim=ylim, ylab = "location",
     xlab = "", main=main, lwd=1, lty=2)
lines(mean.vlo.gaus.fl, type='l', pch=16, ylim=ylim, ylab = "location",
     xlab = "", main=main, col='red', lwd=1, lty=2)
legend("topright", lwd=c(2,2,2,2)/2, lty=c(1,2,1,2), col=cols[c(1,1,2,2)], legend = c("BS (Laplace)","BS (Gaussian)", "FL (Laplace)", "FL (Gaussian)"))
graphics.off()



makejpg(outputdir, "vlo-lapl-minus-gaus.jpg")
main = "Vlo, Laplace minus Gaussian"
ylim = c(-0.1,0.6)
plot(mean.vlo.lapl.bs - mean.vlo.gaus.bs, type='l', pch=16, ylim=ylim, ylab = "location",
     xlab = "", main=main, lwd=1)
lines(mean.vlo.lapl.fl - mean.vlo.gaus.fl, type='l', pch=16, ylim=ylim, ylab = "location",
     xlab = "", main=main, col='red', lwd=1)
abline(h=0)
graphics.off()





## Does the distribution of v^TY differ for locations, when Gaussian and Lapl
## noises are used?

## Generate Y~N(0,1) -->  isolate tests for location i --> See distribution of v^TY 
## Generate Y~Lapl(1/\sqrt{2}) -->  isolate tests for location i --> See distribution of v^TY

## cols=c(rgb(0.5,0,0,0.5), rgb(0,0,0.5,0.5))
## loc = 1:10

## hist(out.gaus$vlo.bs[which(out.gaus$j.bs==loc)], add=FALSE, col=cols[1] ,freq=FALSE, breaks=10, xlim=c(0,5))
## hist(out.lapl$vlo.bs[which(out.lapl$j.bs==loc)], add=TRUE, col =cols[2] , freq=FALSE, breaks=10)

## plot(density(out.gaus$vlo.bs[which(out.gaus$j.bs==loc)]), col=cols[1], xlim=c(0,5))
## lines(density(out.lapl$vlo.bs[which(out.lapl$j.bs==loc)]), col=cols[2])

## summary(out.gaus$vlo.bs[which(out.gaus$j.bs==loc)])
## summary(out.lapl$vlo.bs[which(out.lapl$j.bs==loc)])
## legend("topright", fill=cols, legend = c("gaus","lapl"))


## hist(out.gaus$vlo.bs[which(out.gaus$j.bs==loc)], add=FALSE, col=cols[1] ,freq=FALSE, breaks=10, xlim=c(0,5))
## hist(out.lapl$vlo.bs[which(out.lapl$j.bs==loc)], add=TRUE, col =cols[2] , freq=FALSE, breaks=10)

## hist(out.gaus$vlo.fl[which(out.gaus$j.fl==loc)], add=FALSE, col=cols[1] ,freq=FALSE, breaks=10, xlim=c(0,5))

## hist(out.lapl$vlo.fl[which(out.lapl$j.fl==loc)], add=FALSE, col =cols[1] , freq=FALSE, breaks=10)
## hist(out.lapl$vlo.bs[which(out.lapl$j.bs==loc)], add=TRUE, col =cols[2] , freq=FALSE, breaks=10)




makejpg(outputdir, "vty-comparison-outside.jpg")
cols = c(1,2)
loc=c(1:10,91:100)
plot(density(out.lapl$vlo.bs[which(out.lapl$j.bs%in%loc)]), col=cols[1], xlim=c(0,6),ylim=c(0,2), lwd=2, main="")
lines(density(out.gaus$vlo.bs[which(out.gaus$j.bs%in%loc)]), col=cols[1], xlim=c(0,6),ylim=c(0,2), lwd=2, main="", lty=2)
title(main="vTy for the two distributions")
lines(density(out.lapl$vlo.fl[which(out.lapl$j.fl%in%loc)]), col=cols[2], lwd=2)
lines(density(out.gaus$vlo.fl[which(out.gaus$j.fl%in%loc)]), col=cols[2], xlim=c(0,6),ylim=c(0,2), lwd=2, main="", lty=2)
legend("topright", lwd=c(2,2,2,2), lty=c(1,2,1,2), col=cols[c(1,1,2,2)], legend = c("BS (Laplace)","BS (Gaussian)", "FL (Laplace)", "FL (Gaussian)"))
graphics.off()




makejpg(outputdir, "vty-comparison-inside.jpg")
loc=11:90##c(1:10,91:100)
plot(density(out.lapl$vlo.bs[which(out.lapl$j.bs%in%loc)]), col=cols[1], xlim=c(0,2),ylim=c(0,5), lwd=2, main="")
lines(density(out.gaus$vlo.bs[which(out.gaus$j.bs%in%loc)]), col=cols[1], xlim=c(0,6),ylim=c(0,2), lwd=2, main="", lty=2)
title(main="vTy for the two distributions")
lines(density(out.lapl$vlo.fl[which(out.lapl$j.fl%in%loc)]), col=cols[2], lwd=2)
lines(density(out.gaus$vlo.fl[which(out.gaus$j.fl%in%loc)]), col=cols[2], xlim=c(0,6),ylim=c(0,2), lwd=2, main="", lty=2)
legend("topright", lwd=c(2,2,2,2), lty=c(1,2,1,2), col=cols[c(1,1,2,2)], legend = c("BS (Laplace)","BS (Gaussian)", "FL (Laplace)", "FL (Gaussian)"))
graphics.off()




## makejpg(outputdir, "vlo.jpg")
ylim = c(-0.03,0.7)
plot(meandiffs.bs, type='o', pch=16, ylim=ylim, ylab = "location", xlab = "Vlo(gaus) - Vlo(fl)", main="Vlo shift from Gaus to FL")
lines(meandiffs.fl, pch=16, type='o', ylim=ylim, ylab = "location", xlab = "Vlo(gaus) - Vlo(fl)", main="Vlo shift from Gaus to FL", col='red')
abline(h=0)
graphics.off()




makejpg(outputdir, "middle-vlo-lapl.jpg",500,600)
cols = c(rgb(0.5,0,0,0.5), rgb(0,0.5,0,0.5))
hist(out.lapl$vlo.fl[which(out.lapl$j.fl==50)], ylim=c(0,200), col=cols[1], xlim=c(0,2),
     main="Vup's in the middle (location 50)")
hist(out.lapl$vlo.bs[which(out.lapl$j.bs==50)],add=TRUE, col=cols[2])
legend("topright", fill=cols, legend=c("fused lasso","binary segmentation"))
graphics.off()


makejpg(outputdir, "edge-vlo-lapl.jpg",500,600)
hist(out.lapl$vlo.fl[which(out.lapl$j.fl==10)], ylim=c(0,200), col=cols[1], xlim=c(0,2),
     main="Vup's in the edges (location 10)")
hist(out.lapl$vlo.bs[which(out.lapl$j.bs==10)],add=TRUE, col=cols[2])
legend("topright", fill=cols, legend=c("fused lasso","binary segmentation"))
graphics.off()


makejpg(outputdir, "middle-vlo-gaus.jpg",500,600)
cols = c(rgb(0.5,0,0,0.5), rgb(0,0.5,0,0.5))
hist(out.gaus$vlo.fl[which(out.gaus$j.fl==50)], ylim=c(0,200), col=cols[1], xlim=c(0,2),
     main="Vup's in the middle (location 50)")
hist(out.gaus$vlo.bs[which(out.gaus$j.bs==50)],add=TRUE, col=cols[2])
legend("topright", fill=cols, legend=c("fused lasso","binary segmentation"))
graphics.off()


makejpg(outputdir, "edge-vlo-gaus.jpg",500,600)
hist(out.gaus$vlo.fl[which(out.gaus$j.fl==10)], ylim=c(0,200), col=cols[1], xlim=c(0,2),
     main="Vup's in the edges (location 10)")
hist(out.gaus$vlo.bs[which(out.gaus$j.bs==10)],add=TRUE, col=cols[2])
legend("topright", fill=cols, legend=c("fused lasso","binary segmentation"))
graphics.off()

par(mfrow=c(1,2))



## For fused lasso, the vlo distribution went up a little bit up
summary(out.gaus$vlo.fl[which(out.gaus$j.fl==10)])
summary(out.lapl$vlo.fl[which(out.lapl$j.fl==10)])

## For binary segmentation, the vlo distribution went every so slightly down
summary(out.gaus$vlo.bs[which(out.gaus$j.bs==10)])
summary(out.lapl$vlo.bs[which(out.lapl$j.bs==10)])



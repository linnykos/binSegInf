## Synopsis: Compare Vup and Vlo of one-step fused lasso and bn.aninary segmentation
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
out.lapl = sim(n=100, nrep=20000, err.fun=lapl, seed=NULL)
out.gaus = sim(n=100, nrep=20000, err.fun=rnorm, seed=NULL)



## Make in and out p-value QQ plot
inside = 21:80
outside = c(1:20, 81:100)
outlist = list(out.gaus, out.lapl)
method="bs"## "bs"
for(ii in 1:2){
    out = outlist[[ii]]
    noise = c("gaussian", "laplace")[ii]
    makejpg(outputdir, paste0("pvals-in-and-out-", noise, "-", method, ".jpg"))
    if(method=="bs"){
    qqunif(list(inside=out$p.bs[out$j.bs %in% inside],
                outside=out$p.bs[out$j.bs %in% outside]),
           cols = c(2,1),
           main = noise)
    }
    if(method=="fl"){
    qqunif(list(inside=out$p.fl[out$j.fl %in% inside],
                outside=out$p.fl[out$j.fl %in% outside]),
           cols = c(2,1),
           main = noise)
    }
    legend("topright", col=c("black", "red"), pch=c(1,1), legend=c("outside", "inside"))
    graphics.off()
}



## Plotting vlo vs vty
makejpg(outputdir, "vlo-vs-vty.jpg", 700, 700)
par(mfrow=c(2,2))
outlist = list(out.gaus, out.lapl)
for(ii in 1:2){
    out = outlist[[ii]]
    method = c("gaussian", "laplace")[ii]
    which.inside = which(out$j.bs%in%inside)
    cols = rep("black",20000)
    cols[which.inside] = "red"
    plot(out$vlo.bs ~ out$vty.bs, ## ylim=c(0,6), xlim=c(0,10), 
         ylim=c(0,6),
         col=cols, cex=.3, main="BS")
    title(sub=method)
    abline(v=c(1,2,3,4,5), h=c(1,2,3,4,5), col='grey80')
    legend("topright", col=c("black", "red"), pch=c(1,1), legend=c("outside", "inside"))

    plot(out$vlo.fl ~ out$vty.fl, ## ylim=c(0,6), xlim=c(0,10), 
         ylim=c(0,6),
         col=cols, cex=.3, main="FL")
    title(sub=method)
    abline(v=c(1,2,3,4,5), h=c(1,2,3,4,5), col='grey80')
    legend("topright", col=c("black", "red"), pch=c(1,1), legend=c("outside", "inside"))
}
graphics.off()



## Overall distribution shift
makejpg(outputdir, "vlo-vs-vty-marginal.jpg", 700, 400)

par(mfrow=c(1,2))
vlo.gaus = out.gaus$vlo.bs
vlo.lapl = out.lapl$vlo.bs
plot(density(vlo.lapl), lwd=2, lty=2)
lines(density(vlo.gaus), lwd=2)

vty.gaus = out.gaus$vty.bs
vty.lapl = out.lapl$vty.bs
plot(density(vty.lapl), lwd=2, lty=2)
lines(density(vty.gaus), lwd=2)

graphics.off()




## Overall
makejpg(outputdir, "vty-minus-vlo.jpg", 900, 500)
par(mfrow=c(1,2))
out=out.gaus
n=100
locs = 1:n 
plot(density((out$vty.fl - out$vlo.fl)[out$j.fl%in%locs]), xlim=c(0,6),
     main="FL vty-vlo, overall", ylim=c(0,8), col=2)
out=out.lapl
a = density((out$vty.fl - out$vlo.fl)[out$j.fl%in%locs])
lines(a$y~a$x,col=1)

out=out.gaus
legend("topright", col=c(1,2), lty=c(1,1), legend=c("gaussian","laplace"))
plot(density((out$vty.bs - out$vlo.bs)[out$j.bs%in%locs]), xlim=c(0,6),
     main="BS vty-vlo, overall ", ylim=c(0,8), col=2)
out=out.lapl
a = density((out$vty.bs - out$vlo.bs)[out$j.bs%in%locs]) 
lines(a$y~a$x, col=1)
legend("topright", col=c(1,2), lty=c(1,1), legend=c("lapl","gaus"))
graphics.off()



## This probably intensifies when you isolate to outside
makejpg(outputdir, "vty-minus-vlo-inside.jpg", 1200, 700)
## makejpg(outputdir, "vty-minus-vlo-outside.jpg", 900, 900)
par(mfrow=c(1,2))

## outside = c(1:10, 91:100)
## outside=1:100
outside=11:90

out = out.lapl
diffs.lapl = out$vty.fl - out$vlo.fl
which.outside = which(out$j.fl%in%outside)
diffs.lapl = diffs.lapl[which.outside]
out = out.gaus
diffs.gaus = out$vty.fl - out$vlo.fl
which.outside = which(out$j.fl%in%outside)
diffs.gaus = diffs.gaus[which.outside]
plot(density(diffs.lapl), xlim = c(0,1), main="", ylim=c(0,8), col='red')
lines(density(diffs.gaus), col='red', lty=2)
    legend("topright", col=c(1,2), lty=c(1,1), legend=c("lapl","gaus"))
title(main="FL, vty-vlo, outside")

out = out.lapl
diffs.lapl = out$vty.bs - out$vlo.bs
which.outside = which(out$j.bs%in%outside)
diffs.lapl = diffs.lapl[which.outside]
out = out.gaus
diffs.gaus = out$vty.bs - out$vlo.bs
which.outside = which(out$j.bs%in%outside)
diffs.gaus = diffs.gaus[which.outside]
plot(density(diffs.lapl), xlim = c(0,1), main="")
lines(density(diffs.gaus), col='black', lty=2)
legend("topright", col=c(1,2), lty=c(1,1), legend=c("lapl","gaus"))
title(main="BS, vty-vlo, outside")
legend("topright", lwd=c(2,2,2,2)/2, lty=c(1,2,1,2), col=cols[c(1,1,2,2)],
       legend = c("BS (Laplace)","BS (Gaussian)", "FL (Laplace)", "FL (Gaussian)"))
graphics.off()




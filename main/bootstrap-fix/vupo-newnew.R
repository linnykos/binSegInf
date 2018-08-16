## Synopsis: YET ANOTHER Renewed comparison of Vlo of one-step fused lasso and
## binary segmentation, focusing on different regions
source("helpers.R")
outputdir = "../figures"
la("~/repos/binSegInf/binSegInf")
la("~/repos/genlassoinf/genlassoinf/")

## Run simulation
lapl <- function(n){ rexp(n,rate=sqrt(2)) * sample(c(-1,1),n,replace=TRUE)}
nrep = 20000##500000
out.lapl = sim(n=100, nrep=nrep, err.fun=lapl, seed=NULL)
out.gaus = sim(n=100, nrep=nrep, err.fun=rnorm, seed=NULL)

## Run simulation
loc = c(10)
out = out.gaus

## Binary segmentation
n = 100
a = out$imaxs.bs[which(out$j.bs %in% loc)]
signed.i = table(sapply(a, function(aa)aa[1]))
xtable::xtable(t(signed.i))
makejpg(outputdir, "signed-i.jpg", 800, 500)
par(mfrow=c(1,2))
plot(signed.i, ylim = c(0,45), xlim = c(-n,n), main="Binary segmentation signed i")
abline(v=c(-10,10), col='red', lwd=2)

## Fused lasso
a = out$imaxs.fl[which(out$j.fl %in% loc)]
signed.i = table(sapply(a, function(aa)aa[1]))
xtable::xtable(t(signed.i))
plot(signed.i, ylim = c(0,45), xlim = c(-n,n), main="Fused Lasso signed i")
abline(v=c(-10,10), col='red', lwd=2)
graphics.off()



plot(x=out.gaus$vlo.bs, y=out.gaus$vty.bs)
plot(x=out.gaus$vlo.fl, y=out.gaus$vty.fl)

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

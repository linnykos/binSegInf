## Synopsis: Compare Vlo of one-step fused lasso and binary segmentation
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


out = out.lapl
loc = 1:10

vlo.bs = out$vlo.bs[out$j.bs %in% loc]
vlo.fl = out$vlo.fl[out$j.fl %in% loc]

vty.bs = out$vty.bs[out$j.bs %in% loc]
vty.fl = out$vty.fl[out$j.fl %in% loc]

cols = c(rgb(1/2,0,0,1/2), rgb(0,1/2,0,1/2))

## ylim = c(0, 2)
ylim=c(0,6)
hist(vlo.bs, col=cols[1], freq=FALSE, breaks=20, ylim = ylim)
hist(vlo.fl, add=TRUE, col=cols[2], freq=FALSE, breaks=20)
mean(vlo.bs) - mean(vlo.fl)
legend("topright", fill=cols, legend=c("bs", "fl"))


ylim = c(0, 4)
hist(vty.bs, col=cols[1], freq=FALSE, breaks=20, ylim = ylim)
hist(vty.fl, add=TRUE, col=cols[2], freq=FALSE, breaks=20)
mean(vty.bs) - mean(vty.fl)
legend("topright", fill=cols, legend=c("bs", "fl"))


## I can do a little more; I can get the denominator and the numerator that go
## into Vlo and Vup.
n = 100
mu = rep(0, n)
set.seed(100)
y  = mu + rnorm(n)
out = vlo.surgery(y)

plot(out$vlo.denom.fl/out$vlo.denom.bs)
abline(v=98)
plot(out$vlo.numer.fl/out$vlo.numer.bs)
abline(v=98)

vlo.denoms = sapply(1:nsim, function(isim){
    y  = mu + rnorm(n)
    out = vlo.surgery(y)
    out$vlo.denom.fl })

nsim=1000
outs = lapply(1:nsim, function(isim){
    y  = mu + rnorm(n)
    out = vlo.surgery(y)
    return(out)
})


inside=40:60
outside = 1:10


which.bs.inside = which(sapply(outs, function(out){out$j.bs%in%inside}))
vlo.denom.bs = sapply(outs[which.bs.inside], function(out){out$vlo.denom.bs})

which.fl.inside = which(sapply(outs, function(out){out$j.fl%in%inside}))
vlo.denom.fl =  sapply(outs[which.fl.inside], function(out){out$vlo.denom.fl})

dim(vlo.denom.fl)
dim(vlo.denom.fl)






vlo.denom.bs sapply(outs[which.bs.inside], function(out){out$vlo.denom.bs})
which.bs.outside = which(sapply(outs, function(out){out$j.bs%in%outside}))




out$j.fl
out$j.bs

denom.fl
denom.fl

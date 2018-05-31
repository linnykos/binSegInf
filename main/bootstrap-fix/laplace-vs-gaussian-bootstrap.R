## Synopsis: see non-Gaussian, Laplacian inference as it applies to our scenario.
datadir = "../data"
outputdir = "../output"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))
source(file=file.path("../main/bootstrap-fix/helper.R"))
library(microbenchmark)
library(smoothmest)


## Generate one set of residuals
n = 200
nrepl = 20
mc.cores = 4
m = 50
nsim = 100
lev=1
mn =c(rep(0, m/2), rep(lev, m/2))
pvs.gaus.list = mclapply(1:nrepl, function(irep){
    resid.orig =  rnorm(n)
    pvs.laplace.test = unlist(mclapply(1:nsim, function(isim){
        dosim_bootstrap(resid.orig, size=m, mn=mn, numSteps=10,
                        model = "binseg")
    },mc.cores=mc.cores))
})
pvs.laplace.list = mclapply(1:nrepl, function(irep){
    resid.orig =  rdoublex(n,lambda=1/sqrt(2))
    pvs.laplace.test = unlist(mclapply(1:nsim, function(isim){
        dosim_bootstrap(resid.orig, size=m, mn=mn, numSteps=10,
                        model = "binseg")
    },mc.cores=mc.cores))
})
## filename = "realbootstrap-small-by-step-lev0.Rdata"
filename = "realbootstrap-small-by-step-lev1.Rdata"
save(pvs.gaus.list, pvs.laplace.list, file=file.path(outputdir, filename))


## Plot the results in QQ plots
load(file=file.path(outputdir, filename))
cols = RColorBrewer::brewer.pal(20, "Set4")
plotfilename = makeplotfilename(filename)
jpeg(file=file.path(outputdir, plotfilename), width=500, height=500)
## Four plots
par(mfrow=c(2,2))
qqunif(pvs.gaus.list, type='l', col=cols); title(main="gaussian")
qqunif(unlist(pvs.gaus.list), type='p'); title(main="gaussian, overall")
qqunif(pvs.laplace.list, type='l', col=cols); title(main="laplace")
qqunif(unlist(pvs.laplace.list), type='p'); title(main="laplace, overall")
graphics.off()






## Trying the same with with fused lasso instead of binary segmentation.
n = 200
nrepl = 20
mc.cores = 4
m = 50
nsim = 100
lev=1
mn =c(rep(0, m/2), rep(lev, m/2))
pvs.gaus.list = mclapply(1:nrepl, function(irep){
    resid.orig =  rnorm(n)
    pvs.laplace.test = unlist(mclapply(1:nsim, function(isim){
        dosim_bootstrap(resid.orig, size=m, mn=mn, numSteps=10,
                        model = "fusedlasso")
    },mc.cores=mc.cores))
})
pvs.laplace.list = mclapply(1:nrepl, function(irep){
    resid.orig =  rdoublex(n,lambda=1/sqrt(2))
    pvs.laplace.test = unlist(mclapply(1:nsim, function(isim){
        dosim_bootstrap(resid.orig, size=m, mn=mn, numSteps=10,
                        model = "fusedlasso")
    },mc.cores=mc.cores))
})
filename = "realbootstrap-small-by-step-fusedlasso-lev0.Rdata"
filename = "realbootstrap-small-by-step-fusedlasso-lev1.Rdata"
save(pvs.gaus.list, pvs.laplace.list, file=file.path(outputdir, filename))


## Plot the results in QQ plots
load(file=file.path(outputdir, filename))
cols = RColorBrewer::brewer.pal(20, "Set4")
plotfilename = makeplotfilename(filename)
jpeg(file=file.path(outputdir, plotfilename), width=500, height=500)
## Four plots
par(mfrow=c(2,2))
qqunif(pvs.gaus.list, type='l', col=cols); title(main="gaussian")
qqunif(unlist(pvs.gaus.list), type='p'); title(main="gaussian, overall")
qqunif(pvs.laplace.list, type='l', col=cols); title(main="laplace")
qqunif(unlist(pvs.laplace.list), type='p'); title(main="laplace, overall")
graphics.off()






## Now incorporating ic stopping
n = 200
nrepl = 20
mc.cores = 7
m = 50
nsim = 100
lev=1
mn =c(rep(0, m/2), rep(lev, m/2))
pvs.gaus.list = mclapply(1:nrepl, function(irep){
    resid.orig =  rnorm(n)
    pvs.laplace.test = unlist(mclapply(1:nsim, function(isim){
        dosim_icstop_bootstrap(resid.orig, size=m, mn=mn, model = "binseg")
    },mc.cores=mc.cores))
})
pvs.laplace.list = mclapply(1:nrepl, function(irep){
    resid.orig =  rdoublex(n,lambda=1/sqrt(2))
    pvs.laplace.test = unlist(mclapply(1:nsim, function(isim){
        dosim_icstop_bootstrap(resid.orig, size=m, mn=mn, model = "binseg")
    },mc.cores=mc.cores))
})
filename = "realbootstrap-small-icstop-lev1.Rdata"
## filename = "realbootstrap-small-icstop-lev1.Rdata"
## save(pvs.gaus.list, pvs.laplace.list, file=file.path(outputdir, filename))


## Plot the results in QQ plots
load(file=file.path(outputdir, filename))
cols = RColorBrewer::brewer.pal(20, "Set3")
plotfilename = makeplotfilename(filename)
jpeg(file=file.path(outputdir, plotfilename), width=500, height=500)
## Four plots
par(mfrow=c(2,2))
qqunif(pvs.gaus.list, type='l', col=cols); title(main="gaussian")
qqunif(unlist(pvs.gaus.list), type='p'); title(main="gaussian, overall")
qqunif(pvs.laplace.list, type='l', col=cols); title(main="laplace")
qqunif(unlist(pvs.laplace.list), type='p'); title(main="laplace, overall")
graphics.off()








## Run the example
nsim = 5000
n = 200
mc.cores = 7
lev = 0
mn = c(rep(0, n/2), rep(lev, n/2))
set.seed(0)
resid.laplace = rdoublex(n,lambda=1/sqrt(2))
resid.gaus =  rnorm(n,0,1)
mn = c(rep(0, n/2), rep(lev, n/2))
start.time = Sys.time()
nsim=100
pvs.laplace.test = unlist(mclapply(1:nsim, function(isim){
    printprogress(isim,nsim, start.time=start.time)
    dosim_bootstrap(resid=resid.laplace, mn=mn, numSteps=10, model="binseg")}, mc.cores=mc.cores))
qqunif(pvs.laplace.test)
nsim=100
pvs.gaus.test = unlist(mclapply(1:nsim, function(isim){
    printprogress(isim,nsim, start.time=start.time)
    dosim_bootstrap(resid=resid.gaus, mn=mn, numSteps=10, model="binseg")}, mc.cores=mc.cores))
qqunif(pvs.gaus.test)
filename = "bootstrap-laplace-icstop-fusedlasso-lev0.Rdata"
## filename = "bootstrap-laplace-icstop-fusedlasso-lev1.Rdata"
save(pvs.gaus.test, pvs.laplace.test, file=file.path(outputdir, filename))


lev = 0
n = 200
mn = c(rep(0, n/2), rep(lev, n/2))
y.orig = mn + rdoublex(n,lambda=1/sqrt(2))

bootstrap_sample(y.orig)

if(model=="binseg"){
    h = binSeg_fixedSteps(y, numSteps=numSteps)
    } else {
        h = dualpathSvd2(y=y, D=makeDmat(type="tf",ord=0,m=n), maxsteps=numSteps)
    }
    h.poly = polyhedra(h, numSteps=numSteps)
    vlist = make_all_segment_contrasts_from_cp(cp=h$cp, cp.sign=h$cp.sign, n=n)
    retain = which(sapply(vlist, function(v){all.equal(sum(v*mn),0)==TRUE}))
    if(length(retain)==0) return(NULL)
    vlist = vlist[retain]
    pvs = sapply(vlist, function(v){poly.pval2(y=y, v=v, poly=h.poly, sigma=1)$pv})
    return(pvs)
}









## Now, porting over LAR inference package code and using CPD data matrix


### Fig 6: p values in a high-dimensional model, with 2 nonzero signals
library(selectiveInference)
source("funs.R")
source("boot.R")

set.seed(0)
n = p = 50
cor = "pair"
rho = 0 
s = 2
type = 'gaus'
## bstar = c(rep(4*c(-1,1),length=s),rep(0,p-s))
lev = 4
bstar = c(rep(0,p/2-1), lev, rep(0,p/2))
gen.x <- function(n, center=TRUE){
    D = makeDmat(n,ord=0)
    D = rbind(c(1,rep(0,n-1)), D)
    X = round(ginv(D))
    X = t(t(X) - colMeans(X)) ## crude..
    return(X)
}
nsteps = 3
nboot = 50000
## for (r in 1:R) {
for(isim in 1:nsim){

    ## Generate data 
    x = gen.x(n)
    theta = x %*% bstar
    if (type=='gaus') y = theta + rnorm(n)
    else if (type=='laplace') y = theta + rdoublex(n,lambda=1/sqrt(2))

    ## Fit model
    out = lar(x,y,maxsteps=nsteps,intercept=F,norm=F)
    found[i,,r] = out$action
    a.exact = suppressWarnings(larInf(out,k=nsteps,sigma=1,theta=theta))
    p.exact[i,,r] = a.exact$pv
    u.exact[i,,r] = a.exact$pivot
    for (l in 1:nsteps) {
      a.boot = boot.pval(y,a.exact$vmat[l,],a.exact$vlo[l],a.exact$vup[l],
        nboot,mu=sum(a.exact$vmat[l,]*theta),del.fac)
      p.boot[i,l,r] = a.boot$pv
      u.boot[i,l,r] = a.boot$pivot
    }
    a.sig = suppressWarnings(larInf(out,k=nsteps,sigma=sd(y),theta=theta))
    p.sig[i,,r] = a.sig$pv
    u.sig[i,,r] = a.sig$pivot
  }
}

datdir = "data/"
name = "hi"
save(file=paste0(datdir,".RData"),list=ls())

o

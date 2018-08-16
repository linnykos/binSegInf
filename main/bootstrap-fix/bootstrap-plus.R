## Synopsis: Try out what we're calling "bootstrap-plus" or "bootstrapv2" which
## is to substitute the null distribution of $v^TY$ by some bootstrapped
## residuals around an updated mean (previously grand mean)


library(binSegInf)
library(smoothmest)
outputdir = "../output"
la("~/repos/binSegInf/binSegInf")
la("~/repos/genlassoinf/genlassoinf/")

nrep = 5000
p.tg = p.sub.onejump = p.sub.grandmean = list()

onesim = function(n,  errfun, nboot){

    ## Form data
    mn = rep(0,n)
    y = mn + errfun(n)
    
    ## Fit 1-step binseg algorithm and form contrasts
    obj = binSeg_fixedSteps(y, numSteps=2)
    cp1 = cp1.correction = obj$cp[1]
    cp2 =obj$cp[2]

    ## Make adjusted means
    ym = make_pw_mean(y=y, cp=cp1)
    if(cp1 < n/10 | cp1 > (n-n/10)) cp1.correction = c()
    ym.correction = make_pw_mean(y=y, cp=cp1.correction)

    ## Form contrasts and polyhedon
    vlist = make_all_segment_contrasts(obj)
    which.cp2 = which(abs(as.numeric(names(vlist)))==cp2)
    vlist = vlist[which.cp2]
    G = polyhedra(obj)$gamma
    
    ## Compare p-values of resulting segment test
    p.tg = sapply(vlist, function(v){
        poly.pval(y=y, G=G, u=rep(0,nrow(G)), v, sigma=1)$p
    })
    
    p.sub.onejump = sapply(vlist, function(v){
        pval_plugin_wrapper(y, G, v, nboot=nboot, adjustmean=ym)
    })

    p.sub.onejump.correction = sapply(vlist, function(v){
        pval_plugin_wrapper(y, G, v, nboot=nboot, adjustmean=ym.correction)
    })

    p.sub.grandmean = sapply(vlist, function(v){
        pval_plugin_wrapper(y, G, v, nboot=nboot, adjustmean= rep(mean(y),length(y)))
    })

    return(list(p.tg=p.tg, p.sub.onejump=p.sub.onejump, p.sub.grandmean=p.sub.grandmean,
                p.sub.onejump.correction=p.sub.onejump.correction
                ))
}

## Run corrected version for a single nboot, gaussian noise.
nboot = 10000
results = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim)
    onesim(100, rnorm, nboot)
},mc.cores=6)
save(results, file=file.path(outputdir, "bootstrap-plus-gaus-correction.Rdata"))

## Run corrected version for single nboot, laplace noise
nboot = 10000
results = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim)
    onesim(100, lapl, nboot)
},mc.cores=6)
save(results, file=file.path(outputdir, "bootstrap-plus-lapl-correction.Rdata"))


## Run the experiment for Gaussian data over several |nboot|s.
nsim = 5000
nboots = c(5000, 10000, 20000)
results.by.nboot = list()
for(iboot in 1:length(nboots)){
    nboot = nboots[iboot]
    printprogress(nboot, nboots, fill=TRUE)
    results.by.nboot[[iboot]] = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim)
        onesim(100, rnorm, nboot)
    }, mc.cores=6)
}
save(results.by.nboot, file=file.path(outputdir, "bootstrap-plus-gaus-by-nboot.Rdata"))

## Run the experiment for Laplace data over several |nboot|s.
nsim = 5000
nboots = c(5000, 10000, 20000)
results.by.nboot = list()
for(iboot in 1:length(nboots)){
    nboot = nboots[iboot]
    printprogress(nboot, nboots, fill=TRUE)
    results.by.nboot[[iboot]] = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim)
        onesim(100, lapl, nboot)
    }, mc.cores=6)
}
save(results.by.nboot, file=file.path(outputdir, "bootstrap-plus-lapl-by-nboot.Rdata"))

## Plot either result
## load(file=file.path(outputdir, "bootstrap-plus-gaus-correction.Rdata"))
## load(file=file.path(outputdir, "bootstrap-plus-lapl-correction.Rdata"))
## load(file=file.path(outputdir, "bootstrap-plus-lapl-by-nboot.Rdata"))
load(file=file.path(outputdir, "bootstrap-plus-gaus-by-nboot.Rdata"))

## makejpg(outputdir, "bootstrap-plus-gaus-correction.jpg", mar=c(4.5,4.5,2.5,0.5))
## makejpg(outputdir, "bootstrap-plus-lapl-correction.jpg", mar=c(4.5,4.5,2.5,0.5))
## makejpg(outputdir, "bootstrap-plus-lapl.jpg", mar=c(4.5,4.5,2.5,0.5))
makejpg(outputdir, "bootstrap-plus-gaus.jpg", mar=c(4.5,4.5,2.5,0.5))
ii = 3 ## Just plotting nboot=20000
myresults = results.by.nboot[[ii]]
## myresults = results
p.tg = sapply(myresults, function(a) a$p.tg)
p.sub.onejump = sapply(myresults, function(a) a$p.sub.onejump)
p.sub.grandmean = sapply(myresults, function(a) a$p.sub.grandmean)
p.sub.onejump.correction = sapply(myresults, function(a) a$p.sub.onejump.correction)
qqunif(list(p.tg=p.tg, p.sub.onejump=p.sub.onejump, p.sub.grandmean=p.sub.grandmean),
            ## p.sub.onejump.correction=p.sub.onejump.correction),
       cols=1:4)
## title(main=paste0("number of bootstraps=", nboots[ii]))
title(main="Laplace noise, two types of bootstrap TG p-values after two-step binseg")
graphics.off()

## Proportion of zeros
## par(mfrow=c(1,3))
## for(ii in 1:3){
##     myresults = results.by.nboot[[ii]]
##     ## p.tg = sapply(myresults, function(a) a$p.tg)
##     p.sub.onejump = sapply(myresults, function(a) a$p.sub.onejump)
##     p.sub.onejump = p.sub.onejump[!is.na(p.sub.onejump)]
##     ## p.sub.grandmean = sapply(myresults, function(a) a$p.sub.grandmean)
##     propzero = sum(abs(p.sub.onejump)<1E-10)/length(p.sub.onejump)
##     print(propzero)
## }


# ** How bad is in-sample noise estimation cheating anyway?  As a side issue
# (but important in the big picture), when data is Gaussian but sigma is
# unknown, for saturated tests in changepoint models, we need to estimate them
# in-sample. What damage is being done, if any?  Should we be /seriously/
# concerned, or only partially so?

## Run the experiment with Gaussian noise. 
## nsim = 5000
## nboot = 20000
## start.time=Sys.time()
## results.gaus = mclapply(1:nsim, function(isim){
##     printprogress(isim, nsim, start.time=start.time)
##     onesim(n=100, rnorm, nboot)
## }, mc.cores=3)
## save(results.gaus, file=file.path(outputdir, "bootstrap-plus-gaus.Rdata"))
## load(file=file.path(outputdir, "bootstrap-plus-gaus.Rdata"))
## p.tg = sapply(results.gaus, function(a) a$p.tg)
## p.sub.onejump = sapply(results.gaus, function(a) a$p.sub.onejump)
## p.sub.grandmean = sapply(results.gaus, function(a) a$p.sub.grandmean)
## qqunif(list(p.tg=p.tg, p.sub.onejump=p.sub.onejump, p.sub.grandmean=p.sub.grandmean),
##        cols=1:3)


## ## Run the experiment with Laplac3 noise. 
## nsim = 5000
## nboot = 20000
## start.time=Sys.time()
## results.lapl = mclapply(1:nsim, function(isim){
##     printprogress(isim, nsim, start.time=start.time)
##     onesim(n=100, lapl, nboot)
## }, mc.cores=3)
## save(results.gaus, file=file.path(outputdir, "bootstrap-plus-lapl.Rdata"))

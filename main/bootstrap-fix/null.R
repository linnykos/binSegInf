## Synopsis: produce null p-values from heavy tailed data and using bootstrapped
## null $v^TY$ distribution substitution method.
library(binSegInf)
library(smoothmest)
source("../main/bootstrap-fix/helper.R")
outputdir = "../output"
la("~/repos/binSegInf/binSegInf")
la("~/repos/genlassoinf/genlassoinf/")

## Loda residuals
datadir = "~/repos/binSegInf/data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
samp = resid.cleanmn[-(1:200)]
samp = samp/sd(samp)
rresid <- function(n=NULL, samp=samp){

    ## Notice, n is /not/ used here.
    bootstrap.inds = bootstrap_ind(length(samp), size=length(samp))
    return(samp[bootstrap.inds])
}

## Case 1: All-zero mean null cases:
nrep = 100
## mydf = 3
## out.lapl = simnew(n=100, nrep=nrep, err.fun=lapl, seed=NULL, sigma=1)
## out.gaus = simnew(n=100, nrep=nrep, err.fun=rnorm, seed=NULL, sigma=1)
## out.rt = simnew(n=100, nrep=nrep, err.fun=rt3, seed=NULL, sigma=sqrt(mydf/(mydf-2)))
out.real = simboot(n=length(samp), nrep=nrep,
                   err.fun=rresid, seed=NULL, sigma=1, samp=samp)
save(out.lapl, out.gaus, out.rt, out.real,
     file=file.path(outputdir, "bootstrap-null.Rdata"))

## outputdir = "~/repos/binSegInf/output"
load(file=file.path(outputdir, "bootstrap-null.Rdata"))

## ## Make straightforward p-value plots
makejpg(outputdir, paste0("nullpval-flatmean-lapl.jpg"))
qqunif(out.lapl[11:14], cols=1:4)
title(main="laplace")
graphics.off()

makejpg(outputdir, paste0("nullpval-flatmean-t3.jpg"))
qqunif(out.rt[11:14], cols=1:4)
title(main="t (df=3)")
graphics.off()

makejpg(outputdir, paste0("nullpval-flatmean-gaus.jpg"))
qqunif(out.gaus[11:14], cols=1:4)
title(main="Gaussian")
graphics.off()

makejpg(outputdir, paste0("nullpval-flatmean-real.jpg"))
qqunif(out.real[11:14], cols=1:4)
title(main="Real residuals")
graphics.off()


## Case 2: Nonzero one-jump signal null cases:
nrep = 2500
n = 100
nboot = 10000
rt3.scaled <- function(n){ rt3(n=n,scale=TRUE) }
errfuns <- list(lapl=lapl, rnorm=rnorm, t3=rt3.scaled)
## for(ii in 1:length(errfuns)){
for(ii in 3){
    printprogress(names(errfuns)[ii], names(errfuns),fill=TRUE)
    errfun = errfuns[[ii]]
    for(lev in 1:3){
        printprogress(lev, 1:3, "levs", fill=TRUE)
        p.tg = p.plugin = list()
        mn = c(rep(0, n/2), rep(lev, n/2))
        start.time = Sys.time()
        for(irep in 1:nrep){
    
            ## Generate data
            printprogress(irep, nrep, "replicates", start.time=start.time)
            y = mn + errfun(n)
        
            ## Fit two one-step algorithms and form contrast
            obj = binSeg_fixedSteps(y, numSteps=3)
            vlist = make_all_segment_contrasts(obj)
            tol = 1E-10
        
            ## Return only null contrasts
            whichnull = which(sapply(vlist, function(v) abs(sum(v*mn))<tol ))
            vlist = vlist[whichnull]
            if(length(vlist)==0)next
        
        
            ## Compare p-values of resulting segment test
            G = polyhedra(obj)$gamma
            df = 3
            if(ii==3){ sigma = sqrt(df/(df-2)) } else { sigma = 1 }
            p.tg[[irep]] = sapply(vlist, function(v){
                poly.pval(y=y, G=G, u=rep(0,nrow(G)), v, sigma=sigma)$p
            })
            p.plugin[[irep]] = sapply(vlist, function(v)pval_plugin_wrapper(y, G, v, nboot=nboot, sigma=sigma))
        }
        plist = list(tg=unlist(p.tg), bootstrap=unlist(p.plugin))

        ## Save results
        filename = paste0("bootstrap-nonflat-null-" ,names(errfuns)[ii], "-lev-", lev, ".Rdata")
        save(plist, file=file.path(outputdir, filename))
    
    }
    cat(fill=TRUE)
}

for(ii in 1:length(errfuns)){
    for(lev in 1:3){

        ## Load data
        filename = paste0("bootstrap-nonflat-null-" ,names(errfuns)[ii], "-lev-", lev, ".Rdata")
        load(file=file.path(outputdir, filename))

        ## Make plot
        filename = paste0("nullpval-onejump-lev-", lev, "-", names(errfuns)[ii], ".jpg")
        print(filename)
        makejpg(outputdir,filename)
        qqunif(plist, cols=c(1,2))
        title(main=paste0("Null p-values, One-jump data, n=100, lev=", lev,
                          " and ", names(errfuns)[ii], "noise"))
        graphics.off()
    }
}




## Case 3: Why is the lev=1 giving /so many zeroes/? We want to investigate
## whether it gets better with /more bootstrap replicates? (So that, if we took
## /large/ number of bootstrapped residuals, we would have uniformity.
outputdir = "../output"
la("~/repos/binSegInf/binSegInf")
la("~/repos/genlassoinf/genlassoinf/")
lev = 1
nrep = 3000
n = 10
mn = c(rep(0, n/2), rep(lev, n/2))
nboots = c(500, 1000 ,2000, 5000, 10000, 20000, 100000, 200000)
plist = list()
for(iboot in 1:length(nboots)){
    nboot = nboots[iboot]
    printprogress(nboot, nboots, "number of bootstrap replicates", fill=TRUE)
    errfun = rt3
    p.tg = p.plugin = list()
    for(irep in 1:nrep){
    
        ## Generate data
        printprogress(irep, nrep, "replicates")
        y = mn + errfun(n)
    
        ## Fit two one-step algorithms and form contrast
        obj = binSeg_fixedSteps(y, numSteps=3)
        vlist = make_all_segment_contrasts(obj)
        tol = 1E-10
    
        ## Return only null contrasts
        whichnull = which(sapply(vlist, function(v) abs(sum(v*mn))<tol ))
        vlist = vlist[whichnull]
        if(length(vlist)==0)next
    
        ## Compare p-values of resulting segment test
        G = polyhedra(obj)$gamma
        mydf = 3
        sigma = sqrt(mydf/(mydf-2))
        p.tg[[irep]] = sapply(vlist, function(v){
            poly.pval(y=y, G=G, u=rep(0,nrow(G)), v, sigma=sigma)$p
        })
        p.plugin[[irep]] = sapply(vlist, function(v)
            pval_plugin_wrapper(y, G, v, nboot=nboot, sigma=sigma))
    }
    cat(fill=TRUE)
    plist[[iboot]] = list(p.tg=unlist(p.tg), p.plugin=unlist(p.plugin))
}

filename = "nboot.Rdata"
save(plist, file=file.path(outputdir, filename))


## Make a series of plots that show inferential properties
load(file=file.path(outputdir, filename))

source("../apr12/helpers.R")
makepdf(outputdir,"by-nboot-tg.pdf", width=7, height=7, mar = c(4.5,4.5,0.5,0.5))
qqunif(lapply(plist,function(a)a$p.tg), cols = 1:length(nboots), legend=TRUE)
title(main="lev=1 onejump t3 noise n=100, plain TG p-value distributions \n by number of bootstraps")
graphics.off()


## makepdf(outputdir,"by-nboot-plugin.pdf", width=7, height=7, mar = c(4.5,4.5,0.5,0.5))
makejpg(outputdir,"by-nboot-plugin.jpg", width=700, height=700, mar = c(4.5,4.5,2.5,0.5))
plugin.pvals = lapply(plist,function(a)a$p.plugin)
names(plugin.pvals) = nboots
cols = RColorBrewer::brewer.pal(length(nboots), "Set3")
qqunif(plugin.pvals, cols=cols)
length(plugin.pvals)
title(main="lev=1 onejump t3 noise n=100, bootstrap-substitution p-value distributions \n by number of bootstraps")
graphics.off()


## makepdf(outputdir,"by-nboot-propzero.pdf", width=7, height=7, mar = c(4.5,4.5,0.5,0.5))
makejpg(outputdir,"by-nboot-propzero.jpg", width=700, height=700, mar = c(4.5, 4.5, 2, 0.5))
prop.zeros = sapply(plist,function(a){
    pv = a$p.plugin
    pv = pv[!is.na(pv)]
    sum(abs(pv)<1E-5)/length(pv)
})
names(prop.zeros) = nboots
plot(y=prop.zeros, x=nboots, type='o', xlab = "number of bootstrap replications",
     ylab = "Frequency of pv=0", lwd=2)
title(main="lev=1 onejump t3 noise n=100, \n How often are bootstrap-substitution p-values equal to zero")
graphics.off()


## makepdf(outputdir,"by-nboot-propone.pdf", width=7, height=7, mar = c(4.5, 4.5, 2, 0.5))
makejpg(outputdir,"by-nboot-propone.jpg", width=700, height=700, mar = c(4.5, 4.5, 2, 0.5))
prop.one = sapply(plist,function(a){
    pv = a$p.plugin
    pv = pv[!is.na(pv)]
    sum(abs(pv-1)<1E-5)/length(pv)
})
names(prop.one) = nboots
plot(y=prop.one, x=nboots, type='o', xlab = "number of bootstrap replications", lwd=2,
     ylab = "Frequency of pv=1")
title(main="lev=1 onejump t3 noise n=100, \n How often are bootstrap-substitution p-values equal to one")
graphics.off()






## Why would TG inference be conservative with heavy tails?

## ## Trying flat means now
## source("/media/shyun/Bridge/Dropbox/research/binseginf/notes/apr12/helpers.R")
## rt3 <-function(n){ rt(n, df=3) }
## nrep = 3000##500000
## mydf = 3
## out.lapl = simnew(n=100, nrep=nrep, err.fun=lapl, seed=NULL, sigma=1)
## out.gaus = simnew(n=100, nrep=nrep, err.fun=rnorm, seed=NULL, sigma=1)
## out.rt = simnew(n=100, nrep=nrep, err.fun=rt3, seed=NULL, sigma=sqrt(mydf/(mydf-2)))
## save(out.lapl, out.gaus, out.rt, file=file.path(outputdir, "bootstrap-null.Rdata"))
## qqunif(out.lapl)
## objects(out.lapl)
## names(out.lapl)[11:14] = c("plain.bs", "plain.fl", "bootstrap.bs", "bootstrap.fl")
## qqunif(out.lapl[11:14], cols=c("black","red","grey","pink"))
## ## out.gaus = simnew(n=100, nrep=nrep, err.fun=rnorm, seed=NULL)
## names(out.rt)[11:14] = c("plain.bs", "plain.fl", "bootstrap.bs", "bootstrap.fl")
## qqunif(out.rt[11:14], cols=c("black","red","grey","pink"))

## qqunif(out.rt[11:14], cols=c("black","red","grey","pink"), lty=c(1,1,1,1))



## pp.tg = unlist(p.tg)
## pp.plugin = unlist(p.plugin)
## qqunif(list(tg=pp.tg, plugin=pp.plugin), cols=c("black", "red"))

## ## Make plot
## makejpg(outputdir,paste0("nullpval-onejump-lev-", lev,"-",
##                          names(errfuns)[ii], ".jpg"))
## qqunif(list(plugin=unlist(p.plugin), tg=unlist(p.tg)), cols=c(1,2))
## title(main=paste0("Null p-values, One-jump data, n=10, lev=",lev))
## graphics.off()

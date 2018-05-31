## Synopsis: The sole purpose of this file is to produce a single, very large
## nboot simulation to see if null p-values from heavy tailed data and using
## bootstrapped null $v^TY$ distribution substitution method, is valid.

library(binSegInf)
library(smoothmest)
## source("../apr12/helpers.R")
outputdir = "../output"
la("~/repos/binSegInf/binSegInf")
la("~/repos/genlassoinf/genlassoinf/")
lev = 1
nrep = 3000
n = 10
mn = c(rep(0, n/2), rep(lev, n/2))
nboot = c(100000)
plist = list()
errfun = rt3
p.tg = p.plugin = list()
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
    mydf = 3
    sigma = sqrt(mydf/(mydf-2))
    p.tg[[irep]] = sapply(vlist, function(v){
        poly.pval(y=y, G=G, u=rep(0,nrow(G)), v, sigma=sigma)$p
    })

    p.plugin[[irep]] = sapply(vlist, function(v){

        ## ## Form bootmat times contrast vector in advance
        ## bootmat.times.v

        ## ## Plugin p-value
        pval_plugin_wrapper(y, G, v, nboot=nboot, sigma=sigma)
        ## pval_plugin_wrapper(y, G, v, bootmat.times.v=bootmat.times.v,
        ##                     sigma=sigma)
    })
}
filename = "nboot-somewhatlarge.Rdata"
save(p.tg, p.plugin, file=file.path(outputdir, filename))

outputdir = "../output"
pv = unlist(p.plugin)
pv = pv[!is.na(pv)]
sum(abs(pv)<1E-5)/length(pv)
sum(abs(pv-1)<1E-5)/length(pv)


load(file=file.path(outputdir, filename))
qqunif(unlist(p.tg))
qqunif(unlist(p.plugin))

## ## Make a series of plots
## outputdir = "~/repos/binSegInf/output"
## load(file=file.path(outputdir, filename))
## source("../apr12/helpers.r")
## ## makepdf(outputdir,"by-nboot-tg.pdf", width=7, height=7, mar = c(4.5,4.5,0.5,0.5))
## ## qqunif(lapply(plist,function(a)a$p.tg), cols = 1:length(nboots), legend=TRUE)
## ## title(main="lev=1 onejump t3 noise n=100, plain TG p-value distributions \n by number of bootstraps")
## ## graphics.off()

## ## makepdf(outputdir,"by-nboot-plugin.pdf", width=7, height=7, mar = c(4.5,4.5,0.5,0.5))
## makejpg(outputdir,"by-nboot-plugin.jpg", width=700, height=700, mar = c(4.5,4.5,2.5,0.5))
## plugin.pvals = lapply(plist,function(a)a$p.plugin)
## names(plugin.pvals) = nboots
## qqunif(plugin.pvals, cols = 1:length(nboots))
## title(main="lev=1 onejump t3 noise n=100, bootstrap-substitution p-value distributions \n by number of bootstraps")
## graphics.off()


## ## makepdf(outputdir,"by-nboot-propzero.pdf", width=7, height=7, mar = c(4.5,4.5,0.5,0.5))
## makejpg(outputdir,"by-nboot-propzero.jpg", width=700, height=700, mar = c(4.5, 4.5, 2, 0.5))
## prop.zeros = sapply(plist,function(a){
##     pv = a$p.plugin
##     pv = pv[!is.na(pv)]
##     sum(abs(pv)<1E-5)/length(pv)
## })
## names(prop.zeros) = nboots
## plot(y=prop.zeros, x=nboots, type='o', xlab = "number of bootstrap replications",
##      ylab = "Frequency of pv=0", lwd=2)
## title(main="lev=1 onejump t3 noise n=100, \n How often are bootstrap-substitution p-values equal to zero")
## graphics.off()


## ## makepdf(outputdir,"by-nboot-propone.pdf", width=7, height=7, mar = c(4.5, 4.5, 2, 0.5))
## makejpg(outputdir,"by-nboot-propone.jpg", width=700, height=700, mar = c(4.5, 4.5, 2, 0.5))
## prop.one = sapply(plist,function(a){
##     pv = a$p.plugin
##     pv = pv[!is.na(pv)]
##     sum(abs(pv-1)<1E-5)/length(pv)
## })
## names(prop.one) = nboots
## plot(y=prop.one, x=nboots, type='o', xlab = "number of bootstrap replications", lwd=2,
##      ylab = "Frequency of pv=1")
## title(main="lev=1 onejump t3 noise n=100, \n How often are bootstrap-substitution p-values equal to one")
## graphics.off()





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

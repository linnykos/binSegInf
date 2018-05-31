## Synopsis: Compare real residuals and t3 and laplace distributions in
## quantiles (qqplots), compared to GAUSSIAN.

datadir = "~/repos/binSegInf/data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))

## Make QQ plots to compare (scaled) t3 samples to Gaussian distribution
n = 50000
type = "t3"## type="lapl" ## type="realresid"
for(type in c("t3", "lapl", "realresid")){

    set.seed(0)
    ## Generate samples
    if(type=="t3"){
        samp = rt(n, df=3)
        t3.sd = sqrt(df/(df-2))
        samp = samp/t3.sd
    } else if (type=="lapl"){
        samp = lapl(n)
    } else if (type=="realresid"){
        samp = resid.cleanmn[-(1:200)]
        samp = samp/sd(samp)
    }

    ## Make qqplot
    makejpg(outputdir, paste0("qq-", type,".jpg"), mar=c(4.5,4.5,2.5,0.5))
    qqnorm(samp, main="", ylim = c(-31,31))
    qqline(samp)
    title(main=paste0("Scaled ", type, " samples against N(0,1)"))
    graphics.off()
    
    ## Make histogram
    makejpg(outputdir, paste0("hist-", type,".jpg"), mar=c(4.5,4.5,2.5,0.5))
    hist(samp, breaks=100, main="")
    abline(v= c(-2,-1,0,1,2), col='grey70')
    title(main=paste0("Scaled ", type, " samples and +-{0,1,2}"))
    graphics.off()
}

## Calculate bootstrap-substitute p-values from bootstrapping real residuals
## after $k$ steps, k=2,3,4,5
library(binSegInf)
library(smoothmest)
outputdir = "../output"
la("~/repos/binSegInf/binSegInf")
la("~/repos/genlassoinf/genlassoinf/")

nrep = 5000
p.tg = p.sub.onejump = p.sub.grandmean = list()

onesim = function(n,  errfun, nboot, maxSteps=5){

    ## Form data
    mn = rep(0,n)
    y = mn + errfun(n)
    p.sub.onejump = p.sub.onejump = p.sub.grandmean = rep(NA,maxSteps-1)
    
    for(k in 2:maxSteps){
        ## Fit 1-step binseg algorithm and form contrasts
        obj = binSeg_fixedSteps(y, numSteps=k)
        cp1 =obj$cp[1:(k-1)]
        cp2 =obj$cp[k]
        ym = make_pw_mean(y=y, cp=cp1)

        ## Additionally, we need to make sure that short regions are completely
        ## eliminated... This will take a few minutes to run.

        vlist = make_all_segment_contrasts(obj)
        which.cp2 = which(abs(as.numeric(names(vlist)))==cp2)
        vlist = vlist[which.cp2]
        G = polyhedra(obj)$gamma
        
        ## Compare p-values of resulting segment test
        p.tg[k-1] = sapply(vlist, function(v){
            poly.pval(y=y, G=G, u=rep(0,nrow(G)), v, sigma=1)$p  })
        
        p.sub.onejump[k-1] = sapply(vlist, function(v){
            pval_plugin_wrapper(y, G, v, nboot=nboot, adjustmean=ym)  })

        p.sub.grandmean[k-1] = sapply(vlist, function(v){
            pval_plugin_wrapper(y, G, v, nboot=nboot, adjustmean= rep(mean(y),length(y)))   })
    }

    return(list(p.tg=p.tg, p.sub.onejump=p.sub.onejump, p.sub.grandmean=p.sub.grandmean))
}







## ## Real residuals
## datadir = "~/repos/binSegInf/data"
## filename = "coriell05296.Rdata"
## load(file=file.path(datadir,filename))
## resids.orig = resid.cleanmn[-(1:200)]
## qqnorm(resids.orig)
## qqline(resids.orig)
## ## hist(resids.orig, breaks=100)


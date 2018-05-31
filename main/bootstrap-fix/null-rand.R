## Synopsis: produce *RANDOMIZED* null p-values for using the bootstrap plugin
## of 5.2 in asympinf paper (but with our /known/ residuals instead of $Y-\bar
## Y$)
library(smoothmest)
## source("../apr12/helpers.R")
outputdir = "../output"
la("~/repos/binSegInf/binSegInf")
la("~/repos/genlassoinf/genlassoinf/")

##' Simulation driver for RBS. Takes \code{y.orig} 
onesim_rbs <- function(y=NULL,mn=NULL,bits=1000, fac=1, verbose=FALSE, reduced=FALSE){

    start.time = Sys.time()
    out = inference_bsFs(y=y, max.numSteps=8, consec=2,
                         sigma=sigma, postprocess= TRUE,
                         locs=1:length(y), numIS=10,
                         min.num.things=30,
                         inference.type="pre-multiply",
                         bits=bits, sigma.add=sigma.add,
                         verbose=verbose, start.time=start.time,
                         mn=mn)
    return(out)
}

nsim = 2000
mc.cores = 7
randomized.pvals = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim)

    ## Generate data
    n = 10
    mn = rep(0,n) 
    y = mn + rnorm(n)
    sigma=1
    sigma.add=.2

    ## Fit TG
    a = onesim_rbs(y=y,mn=mn)

    ## Get plugin pvalue if applicable.
    if(is.na(a)){
        return(NULL)
    } else {
        pv = get.plugin.pval(a)
        return(pv)
    }
}, mc.cores=mc.cores)

outputdir = "../output"
save(randomized.pvals, file=file.path(outputdir, "null-rand.Rdata"))



## Plot
outputdir="~/repos/binSegInf/output"
load(file=file.path(outputdir, "null-rand.Rdata"))
makejpg(outputdir, paste0("null-rand.jpg"))
qqunif(unlist(randomized.pvals))
title(main="Gaussian noise, flat signal, n=10, randomized 1-step plugin p-values") 
graphics.off()



## See if the two distributions (bootstrapped, true) are similar?

n = 10
mn = rep(0,n) 
sigma=1
sigma.add=.2


nsim=1000
vy = rep(NA,nsim)
for(isim in 1:nsim){
    y = mn + rnorm(n)
    vy[isim] = v%*%y
}

vystar = rep(NA,nsim)
y = mn + rnorm(n)
for(isim in 1:nsim){
    bootstrap(y)
    v%*%y
}

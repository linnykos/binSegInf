## Synopsis: Using the bootstrap plugin of 5.2 in asympinf paper (but with our
## /known/ residuals instead of $Y-\bar Y$); producing *NONNULL* p-values

library(smoothmest)
source("../main/bootstrap-fix/helpers.R")
outputdir = "../output"
la("~/repos/genlassoinf/genlassoinf/")

## Simulation settings
nrep = 5000
n = 100
nboot = 20000

## Nonzero one-jump signal NONNULL cases:
errfuns <- list(lapl=lapl, rnorm=rnorm, rt=rt)
for(ii in 1:3){
    printprogress(names(errfuns)[ii], names(errfuns),fill=TRUE)
    errfun = errfuns[[ii]]
    for(lev in 1:3){
        printprogress(lev,1:3, "levs")
        p.tg = p.plugin = list()
        mn = c(rep(0, n/2), rep(lev, n/2))
        for(irep in 1:nrep){
    
            ## Generate data
            printprogress(irep, nrep, "replicates")
            y = mn + errfun(n)
        
            ## Fit two one-step algorithms and form contrast
            obj = binSeg_fixedSteps(y, numSteps=3)
            vlist = make_all_segment_contrasts(obj)
            tol = 1E-10
        
            ## Return only null contrasts
            whichnonnull = which(sapply(vlist, function(v) abs(sum(v*mn))>tol ))
            vlist = vlist[whichnonnull]
            if(length(vlist)==0)next
        
            ## Compare polyhedra
            G = polyhedra(obj)$gamma
        
            ## Compare p-values of resulting segment test
            p.tg[[irep]] = sapply(vlist, function(v){
                poly.pval(y=y, G=G, u=rep(0,nrow(G)), v, sigma=1)$p
            })
            p.plugin[[irep]] = sapply(vlist, function(v)pval_plugin_wrapper(y, G, v, nboot=nboot))
        }
    
    }
    cat(fill=TRUE)
}

## ## Make plot
## makejpg(outputdir,paste0("nonnullpval-onejump-lev-", lev,"-",
##                          names(errfuns)[ii], ".jpg"))
## qqunif(list(plugin=unlist(p.plugin), tg=unlist(p.tg)), cols=c(1,2))
## title(main=paste0("Non-null p-values, One-jump data, n=10, lev=",lev))
## graphics.off()

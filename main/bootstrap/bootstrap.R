## Synopsis: run the bootstrap example.

##' Simulation driver for RBS. Takes \code{y.orig} 
onesim_rbs <- function(y.orig, bits=1000, fac=1, verbose=FALSE){

    ## Add bootstrapped residuals around a cleaned mean, with known sigma
    sigma = sd(y.orig[1:200]) * fac
    sigma.add = sigma*0.2
    y = newmn[-(1:200)] + bootstrap_sample(resid.cleanmn[-(1:200)]) * fac
    ## set.seed(1)
    ## lev = 3
    ## n=30
    ## mn = c(rep(0,n/2), rep(lev,n/2))
    ## y = mn + rnorm(n,0,1)
    ## sigma = 1; sigma.add = .2

    start.time = Sys.time()
    out = inference_bsFs(y=y, max.numSteps=15, consec=2,
                         sigma=sigma, postprocess= TRUE,
                         locs=1:length(y), numIS=10,
                         min.num.things=30,
                         inference.type="pre-multiply",
                         bits=bits, sigma.add=sigma.add,
                         verbose=verbose, start.time=start.time,
                         how.close=5,
                         mn=mn)
    return(out)
}


##' Helper to get randomized p-value from simulation results. One caveat.
##' @param myresult is a list that contains \code{y.orig}, \code{vlist}, and
##'     \code{parts.list}. The last thing is a matrix whose column names are
##'     vup/vlo/vty/sigma
##' @return plugin
get.plugin.pval <- function(myresult){

    ## Get bootstrap-centered y's
    y = myresult$y.orig
    y.centered = y - mean(y)
    nboot = 1000
    bootmat = t(sapply(1:nboot, function(iboot){
        y.centered[sample(length(y), size=length(y), replace=TRUE)]
    }))

    ## Produce randomized (importance-weighted) p-values in the plugin bootstrap
    ## distribution
    pvs = rep(NA, length(myresult$vlist))
    for(ii in 1:length(myresult$vlist)){
        parts = myresult$parts.list[[ii]]
        plugin.pvals = apply(parts, 2, function(mycol){
            pval_plugin(Vlo=mycol["vlo"], Vup=mycol["vup"], vty=mycol["vty"],
                        v=myresult$vlist[[ii]], bootmat=bootmat, y=y)
        })
        plugin.weights = apply(parts, 2, function(mycol){
            pval_plugin(Vlo=mycol["vlo"], Vup=mycol["vup"], vty=mycol["vty"],
                        v=myresult$vlist[[ii]], bootmat=bootmat, y=y, weight=TRUE)
        })
        pvs[ii] = sum(plugin.pvals * plugin.weights) / sum(plugin.weights)
    }
    return(pvs)
}

## Synopsis: run the bootstrap example.

##' Simulation driver for RBS. Takes \code{y.orig} 
onesim_rbs <- function(y.orig, bits=1000, fac=1, verbose=FALSE, reduced=FALSE){

    ## Add bootstrapped residuals around a cleaned mean, with known sigma
    sigma = sd(y.orig[1:200]) * fac
    sigma.add = sigma*0.2
    mn = newmn[-(1:200)]
    n = length(mn)
    bootstrap.inds = bootstrap_ind(n, size=n)
    resids.orig = resid.cleanmn[-(1:200)]
    resids = resids.orig[bootstrap.inds]
    y = mn + resids * fac

    start.time = Sys.time()
    out = inference_bsFs(y=y, max.numSteps=15, consec=2,
                         sigma=sigma, postprocess= TRUE,
                         locs=1:length(y), numIS=3,
                         min.num.things=3,
                         inference.type="pre-multiply",
                         bits=bits, sigma.add=sigma.add,
                         verbose=verbose, start.time=start.time,
                         mn=mn,
                         bootstrap.inds=bootstrap.inds)
    return(out)
}


##' Helper to get randomized p-value from simulation results. One caveat.
##' @param myresult is a list that contains \code{y.orig}, \code{vlist}, and
##'     \code{parts.list}. The last thing is a matrix whose column names are
##'     vup/vlo/vty/sigma
##' @return plugin
get.plugin.pval <- function(myresult, nboot=1000){

    ## Get bootstrap-centered y's
    y = myresult$y.orig
    y.centered = y - mean(y)
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

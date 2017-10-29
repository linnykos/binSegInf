##' Helper to get p-values
onejumpsim = function(lev, n, nsim, numSteps, numIS=NULL, randomized){

    ## Basic checks
    if(randomized)assert_that(!is.null(numIS))

    cat("lev=", lev, fill=TRUE)
    sigma = 1

    results = mclapply(1:nsim,function(isim){
        printprogress(isim, nsim)

        ## Generate some data
        mn = c(rep(0,n/2), rep(lev,n/2))
        y = mn + rnorm(n, 0, sigma)

        ## Fit WBS
        numIntervals = n
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
        poly = polyhedra(obj=g$gamma, u=g$u)
        vlist <- make_all_segment_contrasts(g)
        pvs = sapply(vlist, function(v){
            if(randomized){
                return(suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g, sigma=sigma, numIS=numIS)))
            } else {
                return(poly.pval2(y=y, poly=poly, v=v, sigma=sigma)$pv)
            }
        })
        names(pvs) = g$cp*g$cp.sign

        ## Also store other information
        null.true = sapply(vlist, function(v){
            return(v%*%mn == 0)
        })

        return(list(pvs=pvs, null.true=null.true))
    },mc.cores=4)
    cat(fill=TRUE)

    pvs = unlist(lapply(results, function(a)a[["pvs"]]))
    truths = unlist(lapply(results, function(a)a[["null.true"]]))
    return(list(pvs=pvs, truths=truths))
}

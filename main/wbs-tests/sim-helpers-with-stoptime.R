##' Helper to get p-values
dosim_with_stoprule <- function(lev, n, meanfun, nsim, numSteps, numIS=NULL, randomized, mc.cores=4, numIntervals=n,
                                inference.type = "rows", locs=1:n, consec=2){

    ## Basic checks
    if(randomized)assert_that(!is.null(numIS))

    cat("lev=", lev, fill=TRUE)
    sigma = 1

    results = mclapply(1:nsim,function(isim){
        set.seed(isim)
    ## results = lapply(1:nsim,function(isim){
        printprogress(isim, nsim)

        ## Generate some data
        ## mn = c(rep(0,n/2), rep(lev,n/2))
        mn = meanfun(lev,n)
        y = mn + rnorm(n, 0, sigma)
        cumsum.y=cumsum(y)

        ## Fit initial WBS for a generous number of steps
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps,
                                  inference.type='rows')

        ## Get ic-infused polyhedron
        icobj = ic_wrapper(g, sigma=sigma, consec=consec)
        print('stoptime is')
        print(icobj$stoptime)
        if(icobj$stoptime==0) return(NULL) ## If stop time is zero, don't do anything.
        stopped.gamma = do.call(rbind, g$rows.list[1:(icobj$stoptime+consec)])
        stopped.u = rep(0, nrow(stopped.gamma))
        poly = polyhedra(obj=rbind(stopped.gamma, icobj$poly$gamma),
                         u=c(stopped.u, icobj$gamma$u))

        ## Extract changepoints from stopped model and form contrasts
        cp = g$cp[1:icobj$stoptime]
        cp.sign = g$cp.sign[1:icobj$stoptime]
        vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)

        ## Retain only the guys we want
        retain = which((cp %in% locs))
        if(length(retain)==0) return(list(pvs=c(), null.true=c()))


        ## Calculate the p-values
        vlist = vlist[retain] ## Added
        ## pvs = sapply(vlist, function(v){
        pvs = sapply(vlist, function(v){
            if(randomized){
                cumsum.v = cumsum(v)
                return(suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g, sigma=sigma,
                                                        numIS=numIS, inference.type=inference.type,
                                                        cumsum.y=cumsum.y,cumsum.v=cumsum.v,
                                                        stop.time=icobj$stoptime)))
            } else {
                return(poly.pval2(y=y, poly=poly, v=v, sigma=sigma)$pv)
            }
        })
        names(pvs) = (cp*cp.sign)[retain]

        ## Also store other information
        null.true = sapply(vlist, function(v){
            return(v%*%mn == 0)
        })
        ## null.true = null.true[get.rid]

        return(list(pvs=pvs, null.true=null.true))
    },mc.cores=mc.cores)
    ## })
    cat(fill=TRUE)

    pvs = unlist(lapply(results, function(a)a[["pvs"]]))
    truths = unlist(lapply(results, function(a)a[["null.true"]]))
    return(list(pvs=pvs, truths=truths))
}


## Generates one/two-jumped means
onejump <- function(lev,n){c(rep(0,n/2),rep(lev,n/2))}
twojump <- function(lev,n){c(rep(0,n/3),rep(lev,n/3), rep(0,n/3))}
fourjump <- function(lev,n){c(rep(0,n/5), rep(lev,n/5), rep(0,n/5), rep(-2*lev, n/5), rep(0,n/5) )}

##' Helper to get p-values
dosim <- function(lev, n, meanfun, nsim, numSteps, numIS=NULL, randomized, mc.cores=4, numIntervals=n,
                  inference.type = "rows", locs=1:n, better.segment=FALSE, improve.nomass.problem=FALSE){

    ## Basic checks
    if(randomized)assert_that(!is.null(numIS))

    cat("lev=", lev, fill=TRUE)
    sigma = 1

    results = mclapply(1:nsim,function(isim){
        printprogress(isim, nsim)

        ## Generate some data
        mn = meanfun(lev,n)
        y = mn + rnorm(n, 0, sigma)
        cumsum.y=cumsum(y)

        ## Fit WBS
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
        poly = polyhedra(obj=g$gamma, u=g$u)
        if(better.segment){
            vlist <- make_all_segment_contrasts_from_wbs(wbs_obj=g)
        } else {
            vlist <- make_all_segment_contrasts_from_cp(cp=g$cp, cp.sign=g$cp.sign, n=n)
        }

        ## retain only the guys we want
        retain = which((g$cp %in% locs))
        if(length(retain)==0){
            return(list(pvs=c(), null.true=c()))
        }

        ## Get the p-values
        vlist = vlist[retain] ## Added
        pvs = sapply(vlist, function(v){
            if(randomized){
                cumsum.v = cumsum(v)
                return(suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g, sigma=sigma,
                                                        numIS=numIS, inference.type=inference.type,
                                                        cumsum.y=cumsum.y,cumsum.v=cumsum.v,
                                                        improve.nomass.problem=
                                                            improve.nomass.problem)))

            } else {
               return(poly.pval2(y=y, poly=poly, v=v, sigma=sigma)$pv)
            }
        })
        names(pvs) = (g$cp*g$cp.sign)[retain]

        ## Also store other information
        null.true = sapply(vlist, function(v){
            return(v%*%mn == 0)
        })
        ## null.true = null.true[get.rid]

        return(list(pvs=pvs, null.true=null.true))
    },mc.cores=mc.cores)
    cat(fill=TRUE)

    ## results = results[sapply(results,function(a)length(a)==2)]
    pvs = unlist(lapply(results, function(a)a[["pvs"]]))
    truths = unlist(lapply(results, function(a)a[["null.true"]]))
    return(list(pvs=pvs, truths=truths))
}


## Generates one/two-jumped means
onejump <- function(lev,n){c(rep(0,n/2),rep(lev,n/2))}
twojump <- function(lev,n){c(rep(0,n/3),rep(lev,n/3), rep(0,n/3))}
fourjump <- function(lev,n){c(rep(0,n/5), rep(lev,n/5), rep(0,n/5), rep(-2*lev, n/5), rep(0,n/5) )}


dosim_with_stoprule <- function(lev, n, meanfun, nsim, numSteps, numIS=NULL, randomized, mc.cores=4, numIntervals=n,
                                inference.type = "rows", locs=1:n, consec=2, better.segment=FALSE){

    ## Basic checks
    if(randomized)assert_that(!is.null(numIS))

    cat("lev=", lev, fill=TRUE)
    sigma = 1

    results = mclapply(1:nsim,function(isim){
        printprogress(isim, nsim)

        ## Generate some data
        mn = meanfun(lev,n)
        y = mn + rnorm(n, 0, sigma)
        cumsum.y = cumsum(y)

        ## Fit initial WBS for a generous number of steps
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps,
                                  inference.type='rows')

        ## Collect the IC information and polyhedron
        ## Get ic-stopping polyhedron
        ic_obj = get_ic(g$cp, g$y, consec=consec, sigma=sigma, type="bic")
        ## ic_poly = ic_to_poly(ic_obj)
        ic_poly = ic_obj$poly

        ## Check for flag
        if(ic_obj$flag=="normal" ){

            if(!randomized){
                ## Get ic-stopped model selection polyhedron
                stopped.gamma = do.call(rbind, g$rows.list[1:(ic_obj$stoptime+consec)])
                stopped.u = rep(0, nrow(stopped.gamma))
                poly = polyhedra(obj=rbind(stopped.gamma, ic_obj$poly$gamma),
                                 u=c(stopped.u, ic_obj$gamma$u))
            }
        } else {
            return(ic_obj$flag)
        }

        ## Decide on the stoptime
        ## if(icstop){ stoptime = ic_obj$stoptime } else { stoptime = numSteps }
        stoptime  = ic_obj$stoptime

        ## Extract changepoints from stopped model and form contrasts
        cp = g$cp[1:stoptime]
        cp.sign = g$cp.sign[1:stoptime]
        if(better.segment){
            vlist <- make_all_segment_contrasts_from_wbs(wbs_obj=g)
        } else {
            vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
        }

        ## Retain only the changepoints we want results from:
        retain = which((cp %in% locs))
        if(length(retain)==0) return(list(pvs=c(), null.true=c()))

        ## Calculate the p-values
        vlist = vlist[retain]
        pvs = sapply(vlist, function(v){
            if(randomized){
                cumsum.v = cumsum(v)
                return(suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g,
                                                        sigma=sigma,
                                                        numIS=numIS,
                                                        inference.type=inference.type,
                                                        cumsum.y=cumsum.y,
                                                        cumsum.v=cumsum.v,
                                                        stop.time=stoptime+consec,
                                                        ic.poly=ic_poly)))
            } else {
                return(poly.pval2(y=y, poly=poly, v=v, sigma=sigma)$pv)
            }
        })
        names(pvs) = (cp*cp.sign)[retain]

        ## Also store truth
        null.true = sapply(vlist, function(v){
            return(v%*%mn == 0)
        })

        return(list(pvs=pvs, null.true=null.true))
    },mc.cores=mc.cores)
    cat(fill=TRUE)

    ## Classify results
    zero.stop = sum(unlist(sapply(results, function(a)a=="zero.stop")))
    didnt.stop = sum(unlist(sapply(results, function(a)a=="didnt.stop")))
    normal.stop = nsim - didnt.stop - zero.stop
    stop.type = list(zero.stop=zero.stop,didnt.stop=didnt.stop, normal.stop=normal.stop)

    ## Filter for normal results then calculate p-values
    results = results[sapply(results,function(a)length(a)==2)]
    pvs = unlist(lapply(results, function(a)a[["pvs"]]))
    truths = unlist(lapply(results, function(a)a[["null.true"]]))

    return(list(pvs=pvs, truths=truths, stop.type=stop.type))
}



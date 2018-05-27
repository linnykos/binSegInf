##' Does a single randomized wbs (rwbs) inference for a given y
##' @param numIntervals Defaults to \code{length(y)}, or if |intervals| is
##'     provided, the length of that.
##' @param intervals An object of class |intervals|.
##' @return List of simulation results (p-values, locations tested, etc.)
do_rwbs_inference <- function(y=y, max.numSteps=20, numIntervals=length(y),
                              intervals=NULL, consec=2,
                              sigma, postprocess=TRUE, how.close = 5,
                              better.segment=FALSE,
                              locs=1:length(y), numIS=100,
                              inference.type=inference.type,
                              improve.nomass.problem=TRUE, bits=1000,
                              max.numIS=2000,
                              verbose=FALSE, mc.cores=1,
                              min.num.things=30){

    ## Basic checks
    if(!is.null(intervals)){
        if(is.null(intervals)){stop("Provide either |numIntervals| or |intervals|.")}
        numIntervals = intervals$numIntervals
    }

    n=length(y)

    ## Fit initial WBS for a generous number of steps
    if(is.null(intervals) & !is.null(numIntervals)){
        intervals = intervals(numIntervals=numIntervals, n=n)
    }
    g = wildBinSeg_fixedSteps(y, intervals=intervals, numSteps=max.numSteps,
                              inference.type='none')
    cumsum.y = cumsum(y)

    ## Collect the IC information and polyhedron
    ic_obj = get_ic(g$cp, g$y, consec=consec, sigma=sigma, type="bic")
    ic_poly = ic_obj$poly
    stoptime  = ic_obj$stoptime

    if(verbose) cat("stoptime is", stoptime, fill=TRUE)

    ## Check for flag
    if(ic_obj$flag!="normal" ){
        return(NULL)
    }

    ## Extract changepoints from stopped model and declutter
    cp = g$cp[1:stoptime]
    cp.sign = g$cp.sign[1:stoptime]

    if(postprocess){
        cpobj = declutter_new(cp, cp.sign, how.close=how.close)
        cp = abs(cpobj)
        cp.sign = sign(cpobj)
    }

    ## # ## Plot things
    ## plot(y)
    ## abline(v=cp, col='purple',lwd=2)
    ## abline(v=g$cp, col="grey80")
    ## ## abline(v=g$cp[1:ic_obj$stoptime], col='blue', lwd=2)
    ## text(x=g$cp+3, y=rep(1,length(g$cp)), label = 1:length(g$cp))


    ## Form contrasts
    if(better.segment){
        vlist <- make_all_segment_contrasts_from_wbs(wbs_obj=g, cps=cp)
    } else {
        vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=length(y))
    }

    ## Retain only the changepoints we want results from:
    print((abs(as.numeric(names(vlist)))))
    retain = which((abs(as.numeric(names(vlist))) %in% locs))
    if(length(retain)==0) return(list(pvs=c(), null.true=c()))
    vlist = vlist[retain]

    ## Calculate the p-values
    pvs = sapply(1:length(vlist), function(iv){
        printprogress(iv, length(vlist), type = "tests")
        v = vlist[[iv]]
        cumsum.v = cumsum(v)
        pv = suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g,
                                              sigma=sigma,
                                              numIS=numIS,
                                              inference.type=inference.type,
                                              cumsum.y=cumsum.y,
                                              cumsum.v=cumsum.v,
                                              stop.time=stoptime+consec,
                                              ic.poly=ic_poly,
                                              improve.nomass.problem=improve.nomass.problem,
                                              bits=bits,
                                              max.numIS=max.numIS,
                                              mc.cores=mc.cores, min.num.things=min.num.things))
        return(pv)})
    names(pvs) = names(vlist)
    return(list(pvs=pvs, locs.all=cp*cp.sign, locs.retained=as.numeric(names(pvs)), vlist=vlist) )
}



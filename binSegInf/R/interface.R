##' Conducts post-selection inference procedure for additive noise binary
##' segmentation inference -- segment test inference for all detected
##' changepoints -- for a given n-lengthed data vector y, with or without IC
##' stopping.
##' @param added.noise is a manually inputted additive noise vector. This must
##'     be generated from i.i.d. Gaussian noise with \code{sigma} standard
##'     deviation. A rough check is in place, but really just trusting the user
##'     at this point.
##' @param y data vector
##' @param consec The number of rises you'd like as a stopping rule, in terms of
##'     BIC.
##' @param max.numSteps The (maximum) number of algorithm steps you'd like to
##'     do, prior to BIC. If you say \code{icstop=FALSE}, then this is the
##'     actual, fixed number of algorithm steps used.
##' @param icstop Whether to use IC stopping. Defaults to \code{TRUE}.
##' @param verbose Whether to print progress of
##' @param mc.cores Whether to run importance sampler on multiple cores. Not
##'     recommended for normal user.
##' @param start.time Optional argument, for timing things when
##'     \code{verbose=TRUE}
##' @param postprocess \code{TRUE} if you'd like to centroid cluster the
##'     detected changepoints.
##' @param howclose How close would you like points to be in each cluster, for
##'     the centroid clustering, to be.
##' @param return An object of class \code{bsFs} with p-values.
inference_bsFs <- function(y=y, max.numSteps=20, consec=2, sigma, icstop=TRUE,
                                  postprocess=TRUE, locs=1:length(y), numIS=100,
                                  sigma.add = 0.2, bits=50, inference.type=c("rows", "pre-multiply"),
                                  numIntervals=length(y),
                                  max.numIS=2000, verbose=FALSE, min.num.things=10,
                                  added.noise=NULL,
                                  mc.cores=1,
                                  improve.nomass.problem=TRUE,
                                  return.more.things=FALSE,
                                  start.time=NULL,
                                  how.close=5,
                                  whichv = 1){

    ## Basic checks
    inference.type = match.arg(inference.type)
    if(!is.null(added.noise)){
        if( abs(sd(added.noise) - sigma.add) > sigma.add/2){
            stop("Your added noise doesn't match sigma.add well.")
        }
    }

    ## Fit model and get IC information
    n = length(y)
    new.noise = rnorm(length(y),0,sigma.add)
    y.fudged = y + new.noise
    h.fudged = binSeg_fixedSteps(y.fudged, numSteps=max.numSteps)
    h.fudged$y.orig = y
    if(icstop){
        ic_obj = get_ic(h.fudged$cp, h.fudged$y, consec=consec, sigma=sigma+sigma.add, type="bic")
        stoptime = ic_obj$stoptime
        if(ic_obj$flag!="normal"){
            ## return(ic_obj$flag)
            warning(paste0("IC stopping resulted in: ", ic_obj$flag))
            return(h.fudged)
        }
    } else {
        stoptime = max.numSteps
    }

    ## Collect stopped model and postprocess
    cp = h.fudged$cp[1:stoptime]
    cp.sign = h.fudged$cp.sign[1:stoptime]
    if(postprocess){
        cpobj = declutter_new(cp, cp.sign, how.close=how.close)
        cp = abs(cpobj)
        cp.sign = sign(cpobj)
    }

    ## Retain only the changepoints we want results from:
    vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
    retain = which((cp %in% locs))
    if(length(retain)==0) return(list(pvs=c(), null.true=c()))
    vlist = vlist[retain]

    ## Do noise-added inference
    results = lapply(1:length(vlist), function(iv){
        v = vlist[[iv]]
        if(verbose){ printprogress(iv, length(vlist), "segment tests"); cat(fill=TRUE); }
        result = randomize_addnoise(y= y, v=v, sigma=sigma, numIS=numIS,
                                sigma.add=sigma.add, orig.fudged.obj = h.fudged,
                                numSteps = stoptime+consec,
                                ic.poly = ic_obj$poly, bits=bits,
                                inference.type=inference.type,
                                max.numIS=max.numIS, verbose=verbose,
                                mc.cores= mc.cores,
                                improve.nomass.problem=TRUE,
                                return.more.things=TRUE,
                                start.time=start.time,
                                min.num.things=min.num.things
                                )
        return(result)
        if(verbose) cat(fill=TRUE)
    })

    ## names(pvs) = names(vlist)
    ## ## return(list(pvs=pvs, locs.all=cp*cp.sign, locs.retained=as.numeric(names(pvs)), results=results ))
    ## return(results)

    ## The actual end goal of this is to add p-values to the original fudged object
    pvs = sapply(results, function(result)result$pv)
    names(pvs) = names(vlist)
    h.fudged$pvs = pvs
    h.fudged$results = results

    return(h.fudged)
}



##' Does a single randomized wbs (rwbs) inference for a given y
##' @param numIntervals Defaults to \code{length(y)}, or if |intervals| is
##'     provided, the length of that.
##' @param intervals An object of class |intervals|.
##' @param postprocess \code{TRUE} if you'd like to centroid cluster the
##'     detected changepoints.
##' @param howclose How close would you like points to be in each cluster, for
##'     the centroid clustering, to be.
##' @return List of simulation results (p-values, locations tested, etc.)
inference_wbs <- function(y=y, max.numSteps=20, numIntervals=length(y),
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

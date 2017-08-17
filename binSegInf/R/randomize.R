##' Wrapper for doing many wild binary segmentations and computing p-values
##' @param y data vector
##' @param sigma standard deviation of noise
##' @param v Fixed contrast, formed /only/ with the knowledge of the selection
##'     event on \code{y} with some fixed interval, and not from any more
##'     information about \code{y}.
##' @param numIntervals number of WBS intervals you want /each time/. This
##'     should match what you used when applying wild binary segmentation on
##'     your observed dataset.
##' @param numSteps number of steps to take.
##' @param nsim.is Number of importance sampling samples you'd like to
##'     calculate.
##' @param reduce \code{TRUE} if reduced version of polyhedron collecting is to
##'     be used, in polyhedra collecting functions for WBS.
##' @param v Contrast vector.
##' @param augment \code{TRUE} if WBS-FS should be run in augment mode.
##' @example examples/randomized_wildBinSeg_pv-example.R
##' @export
randomized_wildBinSeg_pv <- function(y, sigma, v, cp, numSteps=NULL,
                                     numIntervals, nsim.is, bits=100,
                                     reduce=FALSE, augment=TRUE){

   ## Helper to generate an interval and return /weighted/ inner tg p-value
    get_one <- function(bits=bits){

        .get_cp_from_segment_contrast <- function(v){
            which(dual1d_Dmat(length(v)+2)%*%c(0,v,0)!=0)[2]-1
        }

        .i_covers_cp <- function(i,cp){
            contained = (i$starts <= cp & cp<i$ends)
            return(any(contained))
        }

        ## Generate interval
         ## cp <- .get_cp_from_segment_contrast(v)
        i = generate_intervals(length(y), numIntervals)

        ## Fit new wbs
        obj = wildBinSeg_fixedSteps(y, numSteps, intervals=i, augment=augment)
        ## print(paste('original cp and and new cp are', cp, "and", obj$cp))
        if(length(obj$cp)==0){return(NULL)}

        ## Calculate num & denom of TG
        poly <- polyhedra(obj, reduce=reduce, v=v, sigma=sigma)
        tg = partition_TG(y, poly, v, sigma, nullcontrast=0, bits=bits,reduce=reduce)

        ## Handle case where interval i precludes selection entirely
        ## if(!.i_covers_cp(i,cp)){return(NULL)}

        ## In certain exceptions, we need to add to the denominator but handle the
        ## numerator appropriately.
        tg.behaves.weird = (tg$denom >1 | tg$denom <0 | tg$denom < tg$numer | tg$numer < 0)
        selected.model.was.different = (tg$pv>1 | tg$pv < 0)
        exception <- (!.i_covers_cp(i,cp) | tg.behaves.weird | selected.model.was.different )
        if(exception){
            ## print('exception happened!')
            ## print(c(!.i_covers_cp(i,cp) , tg.behaves.weird , selected.model.was.different))
            if(tg.behaves.weird) print("tg behaves weird!")
            if(selected.model.was.different) print("selected model was different!")
            ## Add this quantity to Wi
            Wi = tg$denom
            ## Add zero to pv
            pv = 0
        } else {
            Wi = tg$denom
            pv = tg$pv
        }


<<<<<<< HEAD
        ## ## Check if tg partition is weird
        ## if(tg$denom >1 | tg$denom <0) print(paste("tg denom was", tg$denom))
        ## tg$denom = min(1, max(tg$denom,0))
        ## if(tg$denom < tg$numer | tg$numer < 0) print(paste("tg numer was", tg$numer))
        ## tg$numer = min(tg$denom, max(tg$numer,0))
=======
        ## Check if tg partition gives any negative or unusual values
        tg$denom = min(1, max(tg$denom,0))
        tg$numer = min(tg$denom, max(tg$numer,0))
>>>>>>> d06b0752ffdfe33b30719a40d72036643ec330ba


        ## return(list(numer = tg$numer, denom = tg$denom, weird=weird))
        return(list(pv=pv, Wi=tg$denom, exception=exception))
    }

    ## Collect weighted p-values and their weights
    pvlist = plyr::rlply(nsim.is, get_one(bit=bits))
    ## pvlist = .filternull(pvlist)

    if(length(pvlist)==0) return(NULL)

    ## Calculate p-value and return
    ## sumNumer = sum(sapply(pvlist, function(nd)nd[["numer"]]))
    ## sumDenom = sum(sapply(pvlist, function(nd)nd[["denom"]]))
    ## pv = sumNumer/sumDenom # sum(unlist(pvmat["numer",]))/ sum(unlist(pvmat["denom",]))

    ## Calculate randomized p-value
    p.vec = sapply(pvlist,function(oneobj)  oneobj[["pv"]] )
    w.vec = sapply(pvlist, function(oneobj) oneobj[["Wi"]])
    exceptions = sapply(pvlist, function(oneobj) oneobj[["exception"]])

    sumNumer = sum(p.vec * w.vec)
    sumDenom = sum(w.vec)
    randomized.pv = sumNumer/sumDenom

    return(list(pv=randomized.pv, p.vec = p.vec, exceptions = exceptions))
}





##' Wrapper for doing many noise-added fused lassos and computing randomized
##' p-value.
##' @param y data vector
##' @param sigma standard deviation of noise
##' @param D penalty matrix.
##' @param v Fixed contrast, formed /only/ with the knowledge of the selection
##'     event on \code{y} with some fixed interval, and not from any more
##'     information about \code{y}.
##' @param numIntervals number of WBS intervals you want /each time/. This
##'     should match what you used when applying wild binary segmentation on
##'     your observed dataset.
##' @param numSteps number of steps to take.
##' @param nsim.is Number of importance sampling samples you'd like to
##'     calculate.
##' @param reduce \code{TRUE} if reduced version of polyhedron collecting is to
##'     be used, in polyhedra collecting functions for WBS.
##' @param v Contrast vector.
##' @param augment \code{TRUE} if WBS-FS should be run in augment mode.
##' @example examples/randomized_wildBinSeg_pv-example.R
##' @export
randomized_genlasso_pv <- function(y, sigma, D, v, numSteps=NULL, numIntervals,
                                   nsim.is, bits=NULL, reduce=FALSE,
                                   augment=TRUE){

    ## Helper to generate an interval and return /weighted/ inner tg p-value
    get_one <- function(bits=bits){

        n = length(y)
        sigmanoise = 0.1 ## sigma*0.1
        noise = rnorm(n,0,1)
        probnoise = prod(sapply(1:length(noise), function(ii) dnorm(noise[ii],0,sigmanoise)))
        ynew = y + noise

        ## Run fused lasso again
        D = genlassoinf::makeDmat(n, type='tf', ord=0)
        fnew = genlassoinf::dualpathSvd2(ynew, D=D, numSteps, approx=TRUE)
        poly <- polyhedra(fnew$Gobj.naive$G, fnew$Gobj.naive$u)
        tg = partition_TG(ynew, poly, v, sigma, nullcontrast=0, bits=bits,reduce=FALSE)

        return(list(numer = tg$numer, denom = tg$denom, probnoise=probnoise ))
    }

    ## Collect weighted p-values and their weights
    pvlist = plyr::rlply(nsim.is, get_one(bit=bits))
    pvlist = .filternull(pvlist)

    if(length(pvlist)==0) return(NULL)

    ## Calculate p-value and return
    ## sumNumer = sum(sapply(pvlist, function(nd)nd[["numer"]]))
    Numers = sapply(pvlist, function(nd)nd[["numer"]])
    Weights = sapply(pvlist, function(nd)nd[["probnoise"]])
    sumNumer = sum(Numers * Weights)

    Denoms = (sapply(pvlist, function(nd)nd[["denom"]]))
    sumDenom = sum(Denoms*Weights)
    ## sumDenom = sum(sapply(pvlist, function(nd)nd[["denom"]]))
    pv = sumNumer/sumDenom # sum(unlist(pvmat["numer",]))/ sum(unlist(pvmat["denom",]))

    return(pv)
}


##' Temporary function to generate |polyhedra| object from fused lasso path
##' output. Mostly exists for formatting purposes.
##' @param obj Output from fused lasso
##' @param reduce If TRUE, then does a Vup/Vlo comparison to see if you should
##'     add (chunks) of rows instead of /all/ of them.
##' @param v a contrast vector, if you want to use smart addition of polyhedra.
##' @param sigma noise level.
##' @param verbose load or not.
##' @return An object of class polyhedra
##' @export
polyhedra_fusedlasso <- function(obj, v=NULL, reduce=FALSE, sigma=NULL,verbose=FALSE,...){

    ## Basic checks
    stopifnot(is_valid.wbsFs(obj))
    if(is.null(v) & reduce) stop("Provide v!")
    if(!is.null(v) & is.null(sigma)) stop("Provide sigma!")

    ## Get all polyhedra
    actual.num.steps = (length(obj$B)-1)

    ## Smartly add rows, if the problem size is big
    if(reduce){
        vup = Inf
        vlo = -Inf
        for(mystep in 1:actual.num.steps){
            newpoly = poly_from_snapshot(obj, mystep, reduce, vup= vup, vlo=vlo,
                                         v=v, sigma=sigma, verbose=verbose)
            vup = newpoly$vup
            vlo = newpoly$vlo
        }
        return(list(v=v,reduce=reduce,obj=obj,vup=vup,vlo=vlo, poly=NULL))

    ## Otherwise, just rbind and add all rows!
    } else {
        all.steps.polys <- lapply(1:actual.num.steps,
                                  function(mystep){
            poly_from_snapshot(obj, mystep, reduce, verbose=verbose)$poly})
        combined.poly = do.call(combine.polyhedra, all.steps.polys)
        return(combined.poly)
    }
}

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
randomized_genlasso_pv <- function(y, sigma, shift, sigma.add, D, v, orig.poly,
                                   numSteps=NULL, numIntervals, nsim.is,
                                   bits=NULL, reduce=FALSE, augment=TRUE){

    ## Helper to generate an interval and return /weighted/ inner tg p-value
    get_one <- function(bits=bits){

        n = length(y)
        new.noise = rnorm(n,0,sigma.add)
        tg = tg_inf(y=y, G=orig.poly$gamma, u=orig.poly$u, shift=new.noise, v=v, sigma=sqrt(sigma^2))
        pv.new = tg$pv
        weight.new = tg$denom

        return
    }

    ## Collect weighted p-values and their weights
    pvlist = plyr::rlply(nsim.is, get_one(bit=bits))
    pvlist = .filternull(pvlist)
    if(length(pvlist)==0) return(NULL)

    ## Calculate p-value and return
    pvs = sapply(pvlist, function(nd)nd[["pv"]])
    denoms = sapply(pvlist, function(nd)nd[["weight"]])
    rtg.pv = sum(pvs*denoms)/sum(denoms)

    return(rtg.pv)
}


##' Function to generate |polyhedra| object from fused lasso path output, from a
##' |path| class object generated from the genlassoinf package.
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




##' Synopsis: randomization wrapper for WBS.
randomize_wbsfs <- function(y, numSteps, numIntervals, numIS = 100, comprehensive=FALSE){

    ## Prepare original set of |numIntervals| intervals
    n = length(y)
    intervals.orig = intervals(numIntervals=numIntervals, n=n, comprehensive=comprehensive)
    if(comprehensive)numIntervals = intervals.orig$numIntervals

    ## Conduct inference once
    y.orig = y
    g.orig = wildBinSeg_fixedSteps(y.orig, numIntervals=numIntervals,
                                   numSteps=numSteps, intervals=intervals.orig)
    poly.orig = polyhedra(obj=g.orig$gamma, u=g.orig$u)
    vlist <- make_all_segment_contrasts(g.orig)
    v = vlist[[1]]

    ##' New PV information based on newly drawn |numIntervals| - |numSteps|
    ##' intervals
    randomized_pv <- function(v, winning.wbs, numIntervals, numSteps){

        ## New intervals added onto old winning intervals
        n = length(v)
        winning_se = rbind(winning.wbs$results[,c("max.s", "max.e")])
        colnames(winning_se) = c("s", "e")
        intervals.new = intervals(numIntervals=numIntervals-numSteps, n=n, existing=winning_se)
        intervals.new = add2(intervals=intervals.new,
                             winning.wbs=g.orig)
        g.new = wildBinSeg_fixedSteps(y=winning.wbs$y,
                    numSteps= numSteps,
                    intervals= intervals.new,
                    mimic=TRUE,
                    wbs.obj=g.orig)
        poly.new = polyhedra(obj=g.new$gamma, u=g.new$u)

        ## Get the conventional p-value, but weigh it by a slightly different
        ## manually constructed weight factor (currently using denom, to see if
        ## it works)
        pvobj = partition_TG(y=y.orig, poly.new, v=v, sigma=sigma)
        pv = pvobj$pv
        weight = pvobj$denom
        ## weight = get_weight()

        info = data.frame(pv=pv,weight=weight)
        return(info)
    }

    if(comprehensive)numIS=1
    pvs = lapply(vlist, function(v){
        things = sapply(1:numIS, function(isim){
            randomized_pv(v=v,
                          winning.wbs=g.orig,
                          numIntervals=numIntervals,
                          numSteps=numSteps
                          )
        })
        pv = sum(unlist(Map('*', things["pv",], things["weight",])))/sum(unlist(things["weight",]))
    })

    names(pvs) = (g.orig$cp) * (g.orig$cp.sign)

    return(pvs)
}


##Iteratively obtain the denominator of the TG statistic
## get_weight <- function(){

##     get_ith_denom <- function(istep, intervals){

##         ## At the i'th step, two things will happen: (1) The i'th win will occur
##         ## in all the eligible intervals at that step, and that the i'th winner
##         ## will have to be smaller than the (i-1)'th winner.

##         ## Get the polyhedra from the (i+1)'th selection
##         ## New intervals added onto old winning intervals
##         winning_se = rbind(winning.wbs$results[,c("max.s", "max.e")])
##         colnames(winning_se) = c("s", "e")
##         intervals.new = intervals(numIntervals=numIntervals-numSteps, n=n, existing=winning_se)
##         intervals.new = add2(intervals=intervals.new,
##                              winning.wbs=g.orig)
##         g.new = wbs(y=winning.wbs$y,
##                     numSteps= numSteps,
##                     intervals= intervals.new,
##                     mimic=TRUE,
##                     wbs.obj=g.orig)
##         poly.new = polyhedra(obj=g.new$rows, u=g.new$u)
##         pv = partition_TG(y=y.orig, poly.new, v=v, sigma=sigma)$pv
##     }
##     }

## }


#######################################
## This is old randomization code #####
#######################################

## ##' Wrapper for doing many wild binary segmentations and computing p-values
## ##' @param y data vector
## ##' @param sigma standard deviation of noise
## ##' @param v Fixed contrast, formed /only/ with the knowledge of the selection
## ##'     event on \code{y} with some fixed interval, and not from any more
## ##'     information about \code{y}.
## ##' @param numIntervals number of WBS intervals you want /each time/. This
## ##'     should match what you used when applying wild binary segmentation on
## ##'     your observed dataset.
## ##' @param numSteps number of steps to take.
## ##' @param nsim.is Number of importance sampling samples you'd like to
## ##'     calculate.
## ##' @param reduce \code{TRUE} if reduced version of polyhedron collecting is to
## ##'     be used, in polyhedra collecting functions for WBS.
## ##' @param v Contrast vector.
## ##' @param augment \code{TRUE} if WBS-FS should be run in augment mode.
## ##' @example examples/randomized_wildBinSeg_pv-example.R
## ##' @export
## randomized_wildBinSeg_pv <- function(y, sigma, v, cp, numSteps=NULL,
##                                      numIntervals, nsim.is, bits=100,
##                                      reduce=FALSE, augment=TRUE){

##    ## Helper to generate an interval and return /weighted/ inner tg p-value
##     get_one <- function(bits=bits){

##         .get_cp_from_segment_contrast <- function(v){
##             which(dual1d_Dmat(length(v)+2)%*%c(0,v,0)!=0)[2]-1
##         }

##         .i_covers_cp <- function(i,cp){
##             contained = (i$starts <= cp & cp<i$ends)
##             return(any(contained))
##         }

##         ## Generate interval; stop when interval i precludes selection
##         ## entirely
##         cp <- .get_cp_from_segment_contrast(v)
##         new.i = generate_intervals(length(y), numIntervals)
##         if(!.i_covers_cp(new.i,cp)){return(NULL)}

##         ## Fit new wbs
##         new.obj = wildBinSeg_fixedSteps(y, numSteps, intervals=new.i, augment=augment)
##         if(length(new.obj$cp)==0){return(NULL)}

##         ## Calculate num & denom of TG
##         new.poly <- polyhedra(new.obj, reduce=reduce, v=v, sigma=sigma)
##         tg = partition_TG(y, new.poly, v, sigma, nullcontrast=0, bits=bits,reduce=FALSE)

##         return(list(pv=tg$pv, Wi=tg$denom))
##     }

##     ## Collect weighted p-values and their weights
##     pvlist = plyr::rlply(nsim.is, get_one(bit=bits))
##     pvlist = .filternull(pvlist)

##     if(length(pvlist)==0)  stop("No inner TG statistics calculated! Try again with bigger |nsim.is|.")

##     ## Calculate randomized p-value without partitioning the num&denom.
##     p.vec = sapply(pvlist, function(oneobj) oneobj[["pv"]] )
##     w.vec = sapply(pvlist, function(oneobj) oneobj[["Wi"]])
##     exceptions = sapply(pvlist, function(oneobj) oneobj[["exception"]])

##     sumNumer = sum(p.vec * w.vec)
##     sumDenom = sum(w.vec)
##     randomized.pv = sumNumer/sumDenom

##     return(list(pv=randomized.pv, p.vec = p.vec, exceptions = exceptions))
## }



## #### Handle exceptions
##         ## Handle three exceptions
##         ## tg.behaves.weird = (tg$denom >1 | tg$denom <0 | tg$denom < tg$numer | tg$numer < 0)
##         ## selected.model.was.different = (tg$pv>1 | tg$pv < 0)
##         ## exception <- (!.i_covers_cp(i,cp) | tg.behaves.weird | selected.model.was.different )
##         ## if(exception){
##         ##     if(tg.behaves.weird) print("tg behaves weird!")
##         ##     if(selected.model.was.different) print("selected model was different!")
##         ##     Wi = tg$denom
##         ##     pv = 0
##         ## } else {
## ## ...
##         ## return(list(numer = tg$numer, denom = tg$denom, weird=weird))
##         ## return(list(pv=pv, Wi=tg$denom, exception=exception))


##         ## ## Check if tg partition gives any negative or unusual values
##         ## tg$denom = min(1, max(tg$denom,0))
##         ## tg$numer = min(tg$denom, max(tg$numer,0))



## ##' Wrapper for doing many wild binary segmentations and computing p-values
## ##' @param y data vector
## ##' @param sigma standard deviation of noise
## ##' @param v Fixed contrast, formed /only/ with the knowledge of the selection
## ##'     event on \code{y} with some fixed interval, and not from any more
## ##'     information about \code{y}.
## ##' @param numIntervals number of WBS intervals you want /each time/. This
## ##'     should match what you used when applying wild binary segmentation on
## ##'     your observed dataset.
## ##' @param numSteps number of steps to take.
## ##' @param nsim.is Number of importance sampling samples you'd like to
## ##'     calculate.
## ##' @param reduce \code{TRUE} if reduced version of polyhedron collecting is to
## ##'     be used, in polyhedra collecting functions for WBS.
## ##' @param v Contrast vector.
## ##' @param augment \code{TRUE} if WBS-FS should be run in augment mode.
## ##' @example examples/randomized_wildBinSeg_pv-example.R
## ##' @export
## randomized_wildBinSeg_pv_new <- function(y, sigma, v, cp, numSteps=NULL,
##                                          numIntervals, winning.s, winning.e,
##                                          nsim.is, bits=100, reduce=FALSE,
##                                          augment=TRUE){

##    ## Helper to generate an interval and return /weighted/ inner tg p-value
##     get_one <- function(bits=bits){
##         .get_cp_from_segment_contrast <- function(v){
##             which(dual1d_Dmat(length(v)+2)%*%c(0,v,0)!=0)[2]-1
##         }

##         .i_covers_cp <- function(i,cp){
##             contained = (i$starts <= cp & cp<i$ends)
##             return(any(contained))
##         }

##         ## Generate (numIntervals-1) intervals, append winning interval
##         new.i.rest = generate_intervals(length(y), numIntervals-1)
##         new.i = add(intervals=new.i.rest,
##                    new.s=winning.s,
##                     new.e=winning.e)

##         ## Fit new wbs
##         new.obj = wildBinSeg_fixedSteps(y, numSteps, intervals=new.i, augment=augment)
##         if(length(new.obj$cp)==0){return(NULL)}

##         ## Calculate num & denom of TG
##         new.poly <- polyhedra(new.obj, reduce=reduce, v=v, sigma=sigma)
##         tg = partition_TG(y, new.poly, v, sigma, nullcontrast=0, bits=bits,reduce=FALSE)

##         return(list(pv=tg$pv, Wi=tg$denom))
##     }

##     ## Collect weighted p-values and their weights
##     pvlist = plyr::rlply(nsim.is, get_one(bit=bits))
##     pvlist = .filternull(pvlist)

##     if(length(pvlist)==0)  stop("No inner TG statistics calculated! Try again with bigger |nsim.is|.")

##     ## Calculate randomized p-value without partitioning the num&denom.
##     p.vec = sapply(pvlist, function(oneobj) oneobj[["pv"]] )
##     w.vec = sapply(pvlist, function(oneobj) oneobj[["Wi"]])
##     exceptions = sapply(pvlist, function(oneobj) oneobj[["exception"]])

##     sumNumer = sum(p.vec * w.vec)
##     sumDenom = sum(w.vec)
##     randomized.pv = sumNumer/sumDenom

##     return(list(pv=randomized.pv, p.vec = p.vec, exceptions = exceptions))
## }

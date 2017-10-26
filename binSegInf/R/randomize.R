##' Function for doing many noise-added fused lassos and computing randomized
##' p-value, from the beginning. It is more desirable to use the wrapper
##' randomized_genlasso(); this function will probably be retired not far from
##' now.
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

##' Wrapper
##' @param pathobj An object from genlassoinf::dualPathSvd2(), of the class
##'     "path".
##' @param sigma.add Standard deviation of noise, for additive noise.
##' @param v contrast vector
##' @param orig.poly original polyhedron to shift.
randomize_genlasso <- function(pathobj, sigma, sigma.add, v, orig.poly,
                                numSteps=NULL, numIntervals, numIS,bits=NULL){

    ## Helper to generate an interval and return /weighted/ inner tg p-value
    get_one <- function(bits=bits){

        n = length(y)
        new.noise = rnorm(n,0,sigma.add)
        tg = partition_TG(y=pathobj$y, poly= polyhedra(obj=orig.poly$gamma,
                                               u=orig.poly$u - orig.poly$gamma%*%new.noise),
                          v=v, sigma=sqrt(sigma^2))
        pv.new = tg$pv
        weight.new = tg$denom
        return(list(pv=pv.new, weight=weight.new))
    }

    ## Collect weighted p-values and their weights
    pvlist = plyr::rlply(numIS, get_one(bit=bits))
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
randomize_wbsfs <- function(v, winning.wbs.obj, numIS = 100, sigma, comprehensive=FALSE){

    numIntervals = winning.wbs.obj$numIntervals
    numSteps = winning.wbs.obj$numSteps

    ##' New PV information based on newly drawn |numIntervals| - |numSteps|
    ##' intervals

    if(comprehensive) numIS=1

    parts = sapply(1:numIS, function(isim){
        rerun_wbs(v=v, winning.wbs.obj=winning.wbs.obj,
                  numIntervals=numIntervals,
                  numSteps=winning.wbs.obj$numSteps,
                  sigma=sigma)
    })
    pv = sum(unlist(Map('*', parts["pv",], parts["weight",])))/sum(unlist(parts["weight",]))

    return(pv)
}

##' Helper for WBSFT randomization, in essence. Rerun WBS to get /new/, single
##' set of denom and numers from a new TG statistic, which is calculated from
##' the /single/ new set of halfspaces that characterize the maximization of the
##' original |winning.wbs.obj| but among /different/ set of
##' |numIntervals|-|numSteps| spaces.
##' @param winning.wbs.obj Original contrast. We call it winning because we will
##'     extract only the winning locations and /those/ winners' enclosing
##'     intervals.
##' @param v test contrast
##' @return A data frame (single row), with "pv" and "weight".
rerun_wbs <- function(winning.wbs.obj, v, numIntervals, numSteps, sigma){

    ## Basic checks
    assert_that(is_valid.wbsFs(winning.wbs.obj))

    ## New intervals added onto old winning intervals
    n = length(v)
    winning_se = rbind(winning.wbs.obj$results[,c("max.s", "max.e")])
    colnames(winning_se) = c("s", "e")
    intervals.new = intervals(numIntervals=numIntervals-numSteps, n=n, existing=winning_se)
    intervals.new = add2(intervals=intervals.new,
                         winning.wbs.obj=winning.wbs.obj)

    ## Create new halfspaces (through |mimic| option)
    g.new = wildBinSeg_fixedSteps(y=winning.wbs.obj$y, numSteps= numSteps,
                                  intervals= intervals.new, mimic=TRUE,
                                  wbs.obj=winning.wbs.obj)
    poly.new = polyhedra(obj=g.new$gamma, u=g.new$u)

    ## Partition TG to denom and numer
    pvobj = partition_TG(y=winning.wbs.obj$y, poly.new, v=v, sigma=sigma, correct.ends=TRUE)
    pv = pvobj$pv
    if(is.nan(pv)) pv=0 ## temporary fix
    weight = pvobj$denom

    ## Special handling so that, if Vup<Vlo, then the weight, which is the prob
    ## along the line trapped in the polyhedron, is zero.
    if(weight<0 | weight>1) weight = 0

    info = data.frame(pv=pv,weight=weight)
    return(info)
}


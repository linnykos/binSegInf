##' Synopsis: noise-added saturated inference, for fused lasso or binary
##' segmentation (really, any method that creates a valid polyhedron and has $cp
##' and $cp.sign)
randomize_addnoise <- function(y, sigma, sigma.add, v, orig.fudged.poly,
                               numSteps=NULL, numIntervals, numIS,bits=50,stopped.poly=NULL, max.numIS = 2000){

    ## New: Get many fudged TG statistics.
    done=FALSE
    pvs = c()
    denoms = c()

    ## Do importance sampling until you have some number of variation..
    while(!done){
        inner.tgs = sapply(1:numIS, function(isim){
            new.noise = rnorm(length(y),0,sigma.add)
            obj.new = partition_TG(y=y, poly=orig.fudged.poly, shift=new.noise,
                                   v=v, sigma=sqrt(sigma^2), bits=bits)
            ## Handle boundary cases
            pv.new = obj.new$pv
            weight.new = obj.new$denom

            ## Handle boundary cases
            if(is.nan(pv.new)) return(c(0,0)) ## Actually not calculable
            if(pv.new>1|pv.new<0)  browser() ## Not sure why this would happen, but anyway!
            if(weight.new<0 | weight.new>1) weight.new=0 ## Nomass problem is to be caught here.
            return(c(pv.new, weight.new))
        })

        rownames(inner.tgs) = c("pv", "denom")
        new.pvs = inner.tgs["pv",]
        new.denoms = inner.tgs["denom",]
        pvs = c(pvs,new.pvs)
        denoms = c(denoms,new.denoms)

        ## increase numIS
        numIS = round(numIS*1.5)

        ## Check if all pvalues are the same, and if so sample more.
        enough.things = any(pvs!=pvs[1])
        reached.limit = numIS > max.numIS
        if(reached.limit | enough.things){ done = TRUE }
    }

    ## Calculate randomized TG statistic
    pv = sum(pvs*denoms)/sum(denoms)
    return(pv)
}


##' Synopsis: randomization wrapper for WBS.
randomize_wbsfs <- function(v, winning.wbs.obj, numIS = 100, sigma,
                            comprehensive=FALSE, inference.type=c("rows", "pre-multiply"),
                            cumsum.y=NULL,cumsum.v=NULL, stop.time=min(winning.wbs.obj$numSteps,
                                                                       length(winning.wbs.obj$cp)),
                            ic.poly=NULL, bits=50, numIS.max=1000,
                            improve.nomass.problem=FALSE, min.num.things=30){

    numIntervals = winning.wbs.obj$numIntervals
    numSteps = winning.wbs.obj$numSteps

    ##' New PV information based on newly drawn |numIntervals| - |numSteps|
    ##' intervals
    if(comprehensive) numIS=1

    done=FALSE
    parts.so.far = cbind(c(Inf,Inf))[,-1,drop=FALSE]
    rownames(parts.so.far) = c("pv", "weight")
    while(!done){
        parts = sapply(1:numIS, function(isim){
            rerun_wbs(v=v, winning.wbs.obj=winning.wbs.obj,
                      numIntervals=numIntervals,
                      numSteps=winning.wbs.obj$numSteps,
                      sigma=sigma,
                      inference.type=inference.type,
                      cumsum.y=cumsum.y,
                      cumsum.v=cumsum.v,
                      stop.time=stop.time,
                      ic.poly=ic.poly)
        })

        parts.so.far = cbind(parts.so.far, parts)
        ## Handling the problem of p-value being NaN/0/1
        things = sum(parts.so.far["weight",]>0)
        enough.things = (things > min.num.things)
        if(!improve.nomass.problem){
            done=TRUE
        }
        if(enough.things){
            done=TRUE
        }
        pv = sum(unlist(Map('*', parts.so.far["pv",], parts.so.far["weight",])))/sum(unlist(parts.so.far["weight",]))
        numIS = numIS*1.5
        if(numIS > numIS.max) done=TRUE
    }
    return(pv)
}

##' Helper for WBSFT randomization, in essence. Rerun WBS to get /new/, singl
##' set of denom and numers from a new TG statistic, which is calculated from
##' the /single/ new set of halfspaces that characterize the maximization of the
##' original |winning.wbs.obj| but among /different/ set of
##' |numIntervals|-|numSteps| spaces.
##' @param winning.wbs.obj Original contrast. We call it winning because we will
##'     extract only the winning locations and /those/ winners' enclosing
##'     intervals.
##' @param v test contrast
##' @return A data frame (single row), with "pv" and "weight".
rerun_wbs <- function(winning.wbs.obj, v, numIntervals, numSteps, sigma,
                      cumsum.y=NULL,cumsum.v=NULL, inference.type, stop.time=numSteps,
                      ic.poly=NULL){

    ## Basic checks
    assert_that(is_valid.wbsFs(winning.wbs.obj))

    ## New intervals added onto old winning intervals
    n = length(v)
    winning_se = rbind(winning.wbs.obj$results[1:stop.time, c("max.s", "max.e")])
    colnames(winning_se) = c("s", "e")
    intervals.new = intervals(numIntervals=numIntervals-stop.time, n=n, existing=winning_se)
    intervals.new = add2(intervals=intervals.new,
                         winning.wbs.obj=winning.wbs.obj,
                         stop.time=stop.time)

    ## Create new halfspaces (through |mimic| option)
    if(inference.type=="rows"){
        g.new = wildBinSeg_fixedSteps(y=winning.wbs.obj$y, numSteps= numSteps,
                                      intervals= intervals.new, mimic=TRUE,
                                      wbs.obj=winning.wbs.obj,
                                      inference.type=inference.type)
        poly.new = polyhedra(obj=g.new$gamma, u=g.new$u)

        ## Partition TG to denom and numer
        pvobj = partition_TG(y=winning.wbs.obj$y, poly.new, v=v, sigma=sigma, correct.ends=TRUE)
        pv = pvobj$pv
        if(is.nan(pv)) pv=0 ## temporary fix
        weight = pvobj$denom


    } else {

        ## new way using new function
        g.new = wildBinSeg_fixedSteps(y=winning.wbs.obj$y, numSteps= numSteps,
                                      intervals= intervals.new, mimic=TRUE,
                                      wbs.obj=winning.wbs.obj,
                                      cumsum.y=cumsum.y,
                                      cumsum.v=cumsum.v,
                                      inference.type="pre-multiply",
                                      stop.time=stop.time,
                                      ic.poly=ic.poly,
                                      v=v)

        ## Calculate TG denom and numer directly
        pvobj = poly_pval_from_inner_products(Gy=g.new$Gy, Gv=g.new$Gv, v=v, y=g.new$y,
                                              sigma=sigma, u=g.new$u, bits=bits)
        pv = pvobj$pv
        if(is.nan(pv)) pv=0 ## temporary fix
        weight = pvobj$denom
    }

    ## Special handling so that, if Vup<Vlo, then the weight, which is the prob
    ## along the line trapped in the polyhedron, is manually assigned zero.
    if(weight<0 | weight>1) weight = 0

    info = data.frame(pv=pv,weight=weight)
    return(info)
}

poly_pval_from_inner_products <- function(Gy,Gv, v,y,sigma,u,bits=50){
    bits=20

    ## Rounding ridiculously small numbers
    Gv[which(abs(Gv)<1E-15)] = 0

    vy = sum(v*y)
    vv = sum(v^2)
    sd = sigma*sqrt(vv)

    rho = Gv / vv
    vec = (u - Gy + rho*vy) / rho
    vlo = suppressWarnings(max(vec[rho>0]))
    vup = suppressWarnings(min(vec[rho<0]))
    vy = max(min(vy, vup),vlo) ##This is the only difference. Should it be here? Yes

    z = Rmpfr::mpfr(vy/sd, precBits=bits)
    a = Rmpfr::mpfr(vlo/sd, precBits=bits)
    b = Rmpfr::mpfr(vup/sd, precBits=bits)

    ## z = vy/sd
    ## a = vlo/sd
    ## b = vup/sd

    if(!(a<=z &  z<=b)){
        warning("F(vlo)<vy<F(vup) was violated, in partition_TG()!")
    }

    ## numer = as.numeric(pnorm(b)-pnorm(z))
    ## denom = as.numeric(pnorm(b)-pnorm(a))
    numer = as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z)))
    denom = as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(a)))
    pv = as.numeric(numer/denom)

    return(list(denom=denom, numer=numer, pv = pv, vlo=vlo, vy=vy, vup=vup))
}

##' Does a single randomized wbs (rwbs) inference for a given y
do_rwbs_inference <- function(y=y, max.numSteps=10, numIntervals=length(y), consec=2,
                             sigma, postprocess=TRUE, better.segment=TRUE,
                             locs=1:length(y), numIS=100,
                             inference.type=inference.type,
                             improve.nomass.problem=TRUE){

    ## Fit initial WBS for a generous number of steps
    g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=max.numSteps,
                              inference.type='none')
    cumsum.y = cumsum(y)

    ## Collect the IC information and polyhedron
    ## Get ic-stopping polyhedron
    ic_obj = get_ic(g$cp, g$y, consec=consec, sigma=sigma, type="bic")
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
    stoptime  = ic_obj$stoptime

    ## Extract changepoints from stopped model and declutter
    cp = g$cp[1:stoptime]
    cp.sign = g$cp.sign[1:stoptime]
    if(postprocess){
        cpobj = declutter(cp=cp, cp.sign)$cp.sign
        cp = cpobj$cp
        cp.sign = cpobj$cp.sign
    }

    ## Form contrasts
    if(better.segment){
        vlist <- make_all_segment_contrasts_from_wbs(wbs_obj=g, cps=cp)
    } else {
        vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
    }

    ## Retain only the changepoints we want results from:
    retain = which((abs(as.numeric(names(vlist))) %in% locs))
    if(length(retain)==0) return(list(pvs=c(), null.true=c()))

    ## Calculate the p-values
    vlist = vlist[retain]
    pvs = sapply(vlist, function(v){
            cumsum.v = cumsum(v)
            return(suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g,
                                                    sigma=sigma,
                                                    numIS=numIS,
                                                    inference.type=inference.type,
                                                    cumsum.y=cumsum.y,
                                                    cumsum.v=cumsum.v,
                                                    stop.time=stoptime+consec,
                                                    ic.poly=ic_poly,
                                                    improve.nomass.problem=improve.nomass.problem)
                                                    ))})
    names(pvs) = (cp*cp.sign)[retain]
    return(list(pvs=pvs, locs=cp[retain]))
}


##' Does a single randomized wbs (rwbs) inference for a given y
do_rfl_inference <- function(y=y, max.numSteps=10, consec=2, sigma,
                             postprocess=TRUE, locs=1:length(y), numIS=100,
                             sigma.add = 0.2, bits=50){


    ## Fit model and get IC information
    n=length(y)
    new.noise = rnorm(length(y),0,sigma.add)
    h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=max.numSteps)
    ic_obj = get_ic(h.fudged$cp, h.fudged$y, consec=consec, sigma=sigma+sigma.add, type="bic")
    stoptime = ic_obj$stoptime
    if(ic_obj$flag!="normal"){return(ic_obj$flag)}

    ## Stopped fudged model (not used because we can't store the polyhedra)
    ## h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=stoptime+consec)
    ## poly.fudged = polyhedra(h.fudged)

    ## Collect stopped model and postprocess
    cp = h.fudged$cp[1:stoptime]
    cp.sign = h.fudged$cp.sign[1:stoptime]
    vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
    ## if(postprocess){
    ##     cpobj = declutter(cp=cp, cp.sign)$cp.sign
    ##     cp = cpobj$cp
    ##     cp.sign = cpobj$cp.sign
    ## }

    ## Retain only the changepoints we want results from:
    retain = which((cp %in% locs))
    if(length(retain)==0) return(list(pvs=c(), null.true=c()))
    vlist = vlist[retain]

    ## Do noise-added inference
    pvs = sapply(vlist, function(v){
        pv = randomize_addnoise(y= y, v=v, sigma=sigma, numIS=numIS,
                                sigma.add=sigma.add, orig.fudged.obj = h.fudged,
                                numSteps = stoptime+consec,
                                ic.poly = ic_obj$poly, bits=bits)
    })

    names(pvs) = (cp*cp.sign)[retain]

    return(list(pvs=pvs, locs=cp[retain]))
}

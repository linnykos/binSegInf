##' Helper to get p-values
dosim <- function(lev, n, meanfun, nsim, numSteps, numIS=NULL, randomized, mc.cores=4, numIntervals=n,
                  inference.type = "rows", locs=1:n, better.segment=FALSE,
                  improve.nomass.problem=FALSE, min.num.things=30, bits=1000){

    ## Basic checks
    if(randomized)assert_that(!is.null(numIS))

    cat("lev=", lev, fill=TRUE)
    sigma = 1

    results = mclapply(1:nsim,function(isim){
        printprogress(isim, nsim)
        set.seed(isim)

        ## Generate some data
        mn = meanfun(lev,n)
        y = mn + rnorm(n, 0, sigma)
        cumsum.y=cumsum(y)

        ## Fit WBS (without polyhedron)
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps,
                                  inference.type="rows")
                                  ## inference.type="fix")
        if(better.segment){
            vlist <- make_all_segment_contrasts_from_wbs(wbs_obj=g)
        } else {
            vlist <- make_all_segment_contrasts_from_cp(cp=g$cp, cp.sign=g$cp.sign, n=n)
        }

        ## Retain only the guys we want
        retain = which((abs(as.numeric(names(vlist))) %in% locs))
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
                                                            improve.nomass.problem, min.num.things=min.num.things,
                                                        bits=bits)))
            } else {
                poly = polyhedra(obj=g$gamma, u=g$u)
               return(poly.pval2(y=y, poly=poly, v=v, sigma=sigma)$pv)
            }
        })
        names(pvs) = names(vlist)

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
                                inference.type = "rows", locs=1:n, consec=2, better.segment=FALSE, improve.nomass.problem=TRUE){

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
                                                        ic.poly=ic_poly,
                                                        improve.nomass.problem=improve.nomass.problem)
                                                        ))
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





## Compare p-values from three methods
dosim_compare <- function(type=c("wbs","fl.nonrand","fl.rand","fl.rand.plus","sbs.rand",
                                 "sbs.nonrand", "wbs.rand", "wbs.nonrand",
                                 "cbs.rand", "cbs.nonrand"),
                          n, lev, numIntervals=n, sigma.add=0.2, numIS=100, meanfun=onejump, visc=NULL, numSteps=1, bits=50){

    type = match.arg(type)
    if(is.null(visc))visc=1:n
    sigma = 1
    mn = meanfun(lev=lev,n=n)
    y = mn + rnorm(n, 0, sigma)
    cumsum.y = cumsum(y)
    inference.type = "pre-multiply"
    improve.nomass.problem = TRUE
    retain=1:n


    if(type == "wbs.rand"){

        ## Fit WBS, test first jump
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
        poly.wbs = polyhedra(obj=g$gamma, u=g$u)
        vlist <- make_all_segment_contrasts(g)

        if(!is.null(visc)){
            retain = which(abs(as.numeric(names(vlist))) %in% visc)
            if(length(retain)==0){
                return(data.frame(pvs=NA, locs=NA))
            }
            vlist = vlist[retain]
        }
        locs = as.numeric(names(vlist))
        retain = which(abs(as.numeric(names(vlist))) %in% visc)

        ## Get the p-values
        pvs = sapply(vlist, function(v){
            cumsum.v = cumsum(v)
            pv = suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g, sigma=sigma,
                                                  numIS=numIS, inference.type=inference.type,
                                                  cumsum.y=cumsum.y,cumsum.v=cumsum.v,
                                                  improve.nomass.problem =improve.nomass.problem,
                                                  bits=bits
                                                  ))
        })
        return(data.frame(pvs=pvs,
                          locs=locs))
    }
    if(type == "wbs.nonrand"){

        ## Fit WBS, test first jump
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
        poly.wbs = polyhedra(obj=g$gamma, u=g$u)
        vlist <- make_all_segment_contrasts(g)
        if(!is.null(visc)){
            retain = which(abs(as.numeric(names(vlist))) %in% visc)
            if(length(retain)==0){
                return(data.frame(pvs=NA, locs=NA))
            }
            vlist = vlist[retain]
        }
        locs = as.numeric(names(vlist))

        pvs = sapply(vlist, function(v){
            cumsum.v = cumsum(v)
            pv = poly.pval2(y=y, poly=poly.wbs, v=v, sigma=sigma,bits=bits)$pv
        })
        return(data.frame(pvs=pvs,
                          locs=locs))
    }

    if(type=="fl.rand"){
        ## Draw new noise
        new.noise = rnorm(n,0,sigma.add)
        ## Fit binseg on fudged data
        D = genlassoinf::makeDmat(n,type='tf',ord=0)
        f.fudged = genlassoinf::dualpathSvd2(y+new.noise, D=D, maxsteps=numSteps, approx=T)
        Gobj.fudged = genlassoinf::getGammat.naive(obj=f.fudged, y=y, condition.step=numSteps) ## Why is this 1?
        poly.fudged = polyhedra(obj=Gobj.fudged$G, u=Gobj.fudged$u)

        ## Get randomized p-value
        vlist <- make_all_segment_contrasts(f.fudged)
        if(!is.null(visc)){
            retain = which(abs(as.numeric(names(vlist))) %in% visc)
            if(length(retain)==0){
                return(data.frame(pvs=NA, locs=NA))
            }
            vlist = vlist[retain]
        }
        locs = as.numeric(names(vlist))

        pvs = sapply(vlist, function(v){
        pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                sigma.add=sigma.add, orig.fudged.poly= poly.fudged, bits=bits)
        })

        return(data.frame(pvs=pvs,
                          locs=locs))
    }

    if(type=="fl.rand.plus"){


        ## Draw new noise
        new.noise = rnorm(n,0,sigma.add)

        ## 2. Get stopped FL p-values
        numSteps = 10
        consec=2
        g.fudged = dualpathSvd2(y + new.noise, maxsteps=numSteps, D = makeDmat(n,ord=0))
        ic_obj = get_ic(g.fudged$cp, g.fudged$y, consec=consec, sigma=sigma, type="bic")
        stoptime = ic_obj$stoptime

        ## Get stopped polyhedron
        if(ic_obj$flag=="normal" ){

            ## Get model selection event polyhedron
            Gobj.fudged = genlassoinf::getGammat.naive(obj=g.fudged, y=y, condition.step=stoptime+consec)
            poly.fudged = polyhedra(obj=Gobj.fudged$G, u=Gobj.fudged$u)

            ## Get ic-stopped model selection polyhedron
            ic_poly = ic_obj$poly

            ## Combine them
            combined.poly = polyhedra(obj = rbind(poly.fudged$gamma, ic_poly$gamma),
                                      u = c(poly.fudged$u, ic_poly$u))

            ## Postprocess and retain vicinity contrasts
            cp = g.fudged$cp[1:stoptime]
            cp.sign = g.fudged$cp.sign[1:stoptime]
            cp = declutter_new(coords=cp, coords.sign=cp.sign, how.close=3)
            cp.sign = sign(cp)
            cp = abs(cp)
            vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
            retain = which(abs(as.numeric(names(vlist))) %in% visc)
            vlist = vlist[retain]

            ## Collect the p-values
            pvs = sapply(vlist, function(v){
                pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                        sigma.add=sigma.add,
                                        orig.fudged.poly= combined.poly)
            })
            locs = as.numeric(names(vlist))
        } else {
            return(data.frame(pvs=NA, locs=NA))
        }
    }

    if(type=="fl.nonrand"){

        ## Get nonrandomized p-value
        D = genlassoinf::makeDmat(n,type='tf',ord=0)
        f.nonfudged = genlassoinf::dualpathSvd2(y, D=D, maxsteps=numSteps, approx=T)
        Gobj.nonfudged = genlassoinf::getGammat.naive(obj=f.nonfudged, y=y, condition.step=numSteps)
        poly.nonfudged = polyhedra(obj=Gobj.nonfudged$G, u=Gobj.nonfudged$u)
        vlist <- make_all_segment_contrasts(f.nonfudged)
        if(!is.null(visc)){
            retain = which(abs(as.numeric(names(vlist))) %in% visc)
            if(length(retain)==0){
                return(data.frame(pvs=NA, locs=NA))
            }
            vlist = vlist[retain]
        }
        locs = as.numeric(names(vlist))

        pvs = sapply(vlist, function(v){
            pv = poly.pval2(y=y, poly=poly.nonfudged, v=v, sigma=sigma, bits=bits)$pv
        })

        return(data.frame(pvs=pvs,
                          locs=locs))
    }




    if(type=="sbs.rand"){

        ## Draw new noise
        new.noise = rnorm(n,0,sigma.add)

        ## Get fudged sbs model
        h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=numSteps)
        poly.fudged = polyhedra(h.fudged)

        ## Get randomized p-value
        vlist <- make_all_segment_contrasts(h.fudged)
        if(!is.null(visc)){
            retain = which(abs(as.numeric(names(vlist))) %in% visc)
            if(length(retain)==0){
                return(data.frame(pvs=NA, locs=NA))
            }
            vlist = vlist[retain]
        }
        locs = as.numeric(names(vlist))
        pvs = sapply(vlist, function(v){
            pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                sigma.add=sigma.add, orig.fudged.poly= poly.fudged, bits=bits)
        })

        return(data.frame(pvs=pvs,
                          locs=locs))
    }


    if(type=="sbs.nonrand"){

        ## Fit binseg on nonfudged data
        h.nonfudged = binSeg_fixedSteps(y, numSteps=numSteps)
        poly.nonfudged = polyhedra(h.nonfudged)

        ## Get randomized p-value
        vlist <- make_all_segment_contrasts(h.nonfudged)
        if(!is.null(visc)){
            retain = which(abs(as.numeric(names(vlist))) %in% visc)
            if(length(retain)==0){
                return(data.frame(pvs=NA, locs=NA))
            }
            vlist = vlist[retain]
        }
        locs = as.numeric(names(vlist))
        pvs = sapply(vlist, function(v){
            pv = poly.pval2(y=y, poly=poly.nonfudged, v=v, sigma=sigma,bits=bits)$pv
        })

        return(data.frame(pvs=pvs,
                          locs=locs))
    }

    if(type=="cbs.rand"){

        ## Draw new noise
        new.noise = rnorm(n,0,sigma.add)

        ## Get fudged sbs model
        stopifnot(numSteps%%2==0)
        h.fudged = circularBinSeg_fixedSteps(y + new.noise, numSteps=numSteps/2)
        poly.fudged = polyhedra(h.fudged)

        ## Get randomized p-value
        vlist <- make_all_segment_contrasts(h.fudged)

        if(!is.null(visc)){
            retain = which(abs(as.numeric(names(vlist))) %in% visc)
            if(length(retain)==0){
                return(data.frame(pvs=NA, locs=NA))
            }
            vlist = vlist[retain]
        }
        locs = as.numeric(names(vlist))

        pvs = sapply(vlist, function(v){
        pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                sigma.add=sigma.add, orig.fudged.poly= poly.fudged,
                                bits=bits)
        })

        return(data.frame(pvs=pvs,
                          locs=locs))
    }

    if(type=="cbs.nonrand"){
        ## Fit cbinseg on nonfudged data
        stopifnot(numSteps%%2==0)
        h.nonfudged = circularBinSeg_fixedSteps(y, numSteps=numSteps/2)
        poly.nonfudged = polyhedra(h.nonfudged)

        ## Get randomized p-value
        vlist <- make_all_segment_contrasts(h.nonfudged)
        if(!is.null(visc)){
            retain = which(abs(as.numeric(names(vlist))) %in% visc)
            if(length(retain)==0){
                return(data.frame(pvs=NA, locs=NA))
            }
            vlist = vlist[retain]
        }
        locs = as.numeric(names(vlist))
        pvs = sapply(vlist, function(v){
        pv = poly.pval2(y=y, poly=poly.nonfudged, v=v, sigma=sigma, bits=bits)$pv
        })

        return(data.frame(pvs=pvs,
                          locs=locs))
    }
}


##' Helper to examine recovery properties of wbs.
dosim_recovery <- function(lev, n, meanfun, nsim, numSteps, numIS=NULL, randomized, mc.cores=4, numIntervals=n,
                  inference.type = "rows", locs=1:n, better.segment=FALSE, improve.nomass.problem=FALSE){

    ## Basic checks
    if(randomized)assert_that(!is.null(numIS))

    sigma=1
    cat("lev=", lev, fill=TRUE)
    results = mclapply(1:nsim,function(isim){
        printprogress(isim, nsim)
        set.seed(isim)

        ## Generate some data
        mn = meanfun(lev,n)
        y = mn + rnorm(n, 0, sigma)

        ## Fit WBS
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)

        ## Calculate proportion
        ## return(g$cp * g$cp.sign)
        return(sum(g$cp %in% locs)/length(g$cp))
    })
    return(mean(unlist(results)))
    ## return(unlist(results))
    ## return(unlist(results))

}


##' Does single fixed-step comparison between fl and bs
dosim_compare_fl_and_bs <- function(){

    ## Simulation settings
    numSteps = 4
    sigma = 1
    meanfun = fourjump
    lev = 2
    mn = meanfun(lev=lev,n=n)
    retain = 1:n
    visc = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
    bits = 1000

    ## Data
    y = mn + rnorm(n, 0, sigma)

    ## Get nonrand binary segmentation
    h.nonfudged = binSeg_fixedSteps(y, numSteps=numSteps)
    poly.nonfudged = polyhedra(h.nonfudged)
    vlist <- make_all_segment_contrasts(h.nonfudged)
    retain = which((h.nonfudged$cp %in% visc))
    if(length(retain)==0){
        pvs.sbs.nonrand = NULL
    } else {
        vlist = vlist[retain]
        pvs.sbs.nonrand = sapply(vlist, function(v){
            pv = poly.pval2(y=y, poly=poly.nonfudged, v=v, sigma=sigma,bits=bits)$pv
        })
    }

    ## Get nonrand fused lasso p-value
    D = genlassoinf::makeDmat(n,type='tf',ord=0)
    f.nonfudged = genlassoinf::dualpathSvd2(y, D=D, maxsteps=numSteps, approx=T)
    Gobj.nonfudged = genlassoinf::getGammat.naive(obj=f.nonfudged, y=y, condition.step=1)
    poly.nonfudged = polyhedra(obj=Gobj.nonfudged$G, u=Gobj.nonfudged$u)
    vlist <- make_all_segment_contrasts(f.nonfudged)
    retain = which((f.nonfudged$cp %in% visc))
    if(length(retain)==0){
        ## return(data.frame(pvs=NA, locs=NA))
        pvs.fl.nonrand = NULL
    } else {
        vlist = vlist[retain]
        pvs.fl.nonrand = sapply(vlist, function(v){
            pv = poly.pval2(y=y, poly=poly.nonfudged, v=v, sigma=sigma, bits=bits)$pv
        })
    }

    return(list(pvs.fl.nonrand=pvs.fl.nonrand,
                pvs.sbs.nonrand=pvs.sbs.nonrand))
}


## Synopsis 2: Compare FL and BS (nonrand) version with IC stopping and with
## decluttering
dosim_compare_fl_and_bs_with_stoprule_and_decluttering <- function(){

    ## Simulation settings
    n = 200  ## n=50
    consec = 2
    meanfun = fourjump
    lev = 1
    mn = meanfun(lev,n)
    sigma = 1
    y = mn + rnorm(n, 0, sigma)
    cumsum.y = cumsum(y)
    sigma.add = 0.2
    numSteps = 15
    numIS = 100
    new.noise = rnorm(n,0,sigma.add)
    visc = sapply(1:4, function(ii){ii/5*n + c(-2,-1,0,1,2)})
    improve.nomass.problem = TRUE
    inference.type="pre-multiply"

    ## 1. Get stopped SBS p-values
    h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=numSteps)
    ic_obj = get_ic(h.fudged$cp, h.fudged$y, consec=consec, sigma=sigma, type="bic")
    stoptime1 = ic_obj$stoptime

    ## Get stopped polyhedron
    if(ic_obj$flag=="normal" ){

        ## Get model selection event polyhedron
        poly.fudged = polyhedra(h.fudged, numSteps = stoptime1+consec)

        ## Get ic-stoppage polyhedron
        ic_poly = ic_obj$poly

        ## Combine them
        combined.poly = polyhedra(obj = rbind(poly.fudged$gamma, ic_poly$gamma),
                                  u = c(poly.fudged$u, ic_poly$u))

        ## Postprocess and retain vicinity contrasts
        cp = h.fudged$cp[1:stoptime1]
        cp.sign = h.fudged$cp.sign[1:stoptime1]
        vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
        numtests.before.retain.sbs <- length(vlist)
        retain = which(abs(as.numeric(names(vlist))) %in% visc)
        vlist = vlist[retain]

        ## Collect the p-values
        pvs.sbs = sapply(vlist, function(v){
            pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                    sigma.add=sigma.add,
                                    orig.fudged.poly=combined.poly)
        })
    } else {
        pvs.sbs=NULL
    }

    ## 2. Get stopped FL p-values
    g.fudged = dualpathSvd2(y + new.noise, maxsteps=numSteps, D = makeDmat(n,ord=0))
    ic_obj = get_ic(g.fudged$cp, g.fudged$y, consec=consec, sigma=sigma, type="bic")
    stoptime2 = ic_obj$stoptime

    ## Get stopped polyhedron
    if(ic_obj$flag=="normal" ){

        ## Get model selection event polyhedron
        Gobj.fudged = genlassoinf::getGammat.naive(obj=g.fudged, y=y, condition.step=stoptime2+consec)
        poly.fudged = polyhedra(obj=Gobj.fudged$G, u=Gobj.fudged$u)

        ## Get ic-stopped model selection polyhedron
        ic_poly = ic_obj$poly

        ## Combine them
        combined.poly = polyhedra(obj = rbind(poly.fudged$gamma, ic_poly$gamma),
                                  u = c(poly.fudged$u, ic_poly$u))

        ## Postprocess and retain vicinity contrasts
        cp = g.fudged$cp[1:stoptime1]
        cp.sign = g.fudged$cp.sign[1:stoptime1]
        cp = declutter_new(coords=cp, coords.sign=cp.sign, how.close=3)
        cp.sign = sign(cp)
        cp = abs(cp)
        vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
        numtests.before.retain.fl <- length(vlist)
        retain = which(abs(as.numeric(names(vlist))) %in% visc)
        vlist = vlist[retain]

        ## Collect the p-values
        pvs.fl = sapply(vlist, function(v){
            pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                    sigma.add=sigma.add,
                                    orig.fudged.poly= combined.poly)
        })
    } else {
        pvs.fl = NULL
   }

    return(list(pvs.sbs=pvs.sbs,
                numtests.before.retain.sbs=numtests.before.retain.sbs,
                pvs.fl=pvs.fl,
                numtests.before.retain.fl=numtests.before.retain.fl
                ))
}



##' Helper function to compute p-values from three methods
dosim_compare <- function(type=c("fl","fl.noisy","fl.noisy.plus",
                                 "sbs.noisy.nonmarg", "sbs.noisy", "sbs",
                                 "wbs.marg", "wbs", "cbs", "cbs.noisy"), n, lev,
                          numIntervals=n, sigma.add=0.2, numIS=100,
                          meanfun=onejump, visc=NULL, numSteps=1, bits=1000,
                          max.numIS=2000, verbose=FALSE, min.num.things=30,
                          sigma = 1){

    type = match.arg(type)
    if(is.null(visc))visc=1:n
    mn = meanfun(lev=lev,n=n)
    y = mn + rnorm(n, 0, sigma)
    cumsum.y = cumsum(y)
    inference.type = "pre-multiply"


    if(type == "wbs.marg"){

        ## Fit WBS, test first jump
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
        poly.wbs = polyhedra(obj=g$gamma, u=g$u)
        vlist <- make_all_segment_contrasts(g)
        vlist <- filter_vlist(vlist, visc)
        locs = as.numeric(names(vlist))
        retain = which(abs(as.numeric(names(vlist))) %in% visc)

        ## Get the p-values
        pvs = sapply(vlist, function(v){
            cumsum.v = cumsum(v)
            pv = suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g, sigma=sigma,
                                                  numIS=numIS, inference.type=inference.type,
                                                  cumsum.y=cumsum.y,cumsum.v=cumsum.v,
                                                  bits=bits, max.numIS=max.numIS, verbose=verbose,
                                                  min.num.things=min.num.things
                                                  ))
        })
        return(data.frame(pvs=pvs,
                          locs=locs))
    }
    if(type == "wbs"){

        ## Fit WBS, test first jump
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
        poly.wbs = polyhedra(obj=g$gamma, u=g$u)
        vlist <- make_all_segment_contrasts(g)
        vlist <- filter_vlist(vlist, visc)
        locs = as.numeric(names(vlist))

        pvs = sapply(vlist, function(v){
            cumsum.v = cumsum(v)
            pv = poly.pval2(y=y, poly=poly.wbs, v=v, sigma=sigma,bits=bits)$pv
        })
        return(data.frame(pvs=pvs,
                          locs=locs))
    }

    if(type=="fl.noisy"){

        ## Draw new noise
        new.noise = rnorm(n,0,sigma.add)

        ## Fit FL on fudged data
        D = genlassoinf::makeDmat(n,type='tf',ord=0)
        f.fudged = genlassoinf::dualpathSvd2(y+new.noise, D=D, maxsteps=numSteps, approx=T)
        Gobj.fudged = genlassoinf::getGammat.naive(obj=f.fudged, y=y, condition.step=numSteps) ## Why is this 1?
        poly.fudged = polyhedra(obj=Gobj.fudged$G, u=Gobj.fudged$u)

        ## Get randomized p-value
        vlist <- make_all_segment_contrasts(f.fudged)
        vlist <- filter_vlist(vlist, visc)
        locs = as.numeric(names(vlist))

        pvs = sapply(vlist, function(v){
            pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                    sigma.add=sigma.add,
                                    orig.fudged.poly= poly.fudged, bits=bits,
                                    max.numIS=max.numIS, verbose=verbose,
                                    min.num.things=min.num.things)$pv})

        return(data.frame(pvs=pvs,
                          locs=locs))
    }

    if(type=="fl.noisy.plus"){


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
            vlist <- filter_vlist(vlist, visc)

            ## Collect the p-values
            pvs = sapply(vlist, function(v){
                pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                        sigma.add=sigma.add,
                                        orig.fudged.poly= combined.poly,
                                        bits=bits, max.numIS=max.numIS,
                                        verbose=verbose,
                                        min.num.things=min.num.things
                                        )$pv})
            locs = as.numeric(names(vlist))
            return(data.frame(pvs=pvs, locs=locs))
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
        vlist <- filter_vlist(vlist, visc)
        locs = as.numeric(names(vlist))

        pvs = sapply(vlist, function(v){
            pv = poly.pval2(y=y, poly=poly.nonfudged, v=v, sigma=sigma, bits=bits)$pv
        })

        return(data.frame(pvs=pvs,
                          locs=locs))
    }



    if(type=="sbs.rand.nonmarg"){

        ## Draw new noise to obfuscate model
        new.noise = rnorm(n,0,sigma.add)

        ## Get fudged sbs model
        h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=numSteps)
        poly.fudged = polyhedra(h.fudged)

        ## Make contrasts
        vlist <- make_all_segment_contrasts(h.fudged)
        vlist <- filter_vlist(vlist, visc)
        locs = as.numeric(names(vlist))

        ## Calculate TG statistic under /shifted/ poly from another drawn noise.
        pvs = sapply(vlist, function(v){
            obj.new = partition_TG(y=y, poly=poly.fudged, shift=new.noise,
                                   v=v, sigma=sqrt(sigma^2), bits=bits)
            pv.new = obj.new$pv
            return(pv.new)
        })

        ## ## Alternative
        ## pvs = sapply(vlist, function(v){
        ##     pv = poly.pval2(y=y, poly=poly.fudged, v=v, sigma=sigma,bits=bits)$pv

        ## })
        return(data.frame(pvs=pvs,
                          locs=locs))
    }


    if(type=="sbs.noisy"){

        ## Draw new noise
        new.noise = rnorm(n,0,sigma.add)

        ## Get fudged sbs model
        h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=numSteps)
        poly.fudged = polyhedra(h.fudged)

        ## Get randomized p-value
        vlist <- make_all_segment_contrasts(h.fudged)
        vlist <- filter_vlist(vlist, visc)
        locs = as.numeric(names(vlist))
        pvs = sapply(vlist, function(v){
            pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                    sigma.add=sigma.add,
                                    orig.fudged.poly=poly.fudged, bits= bits,
                                    max.numIS=max.numIS, verbose=verbose,
                                    min.num.things=min.num.things)$pv})
        return(data.frame(pvs=pvs,
                          locs=locs))
    }


    if(type=="sbs"){

        ## Fit binseg on nonfudged data
        h.nonfudged = binSeg_fixedSteps(y, numSteps=numSteps)
        poly.nonfudged = polyhedra(h.nonfudged)

        ## Get randomized p-value
        vlist <- make_all_segment_contrasts(h.nonfudged)
        vlist <- filter_vlist(vlist, visc)
        locs = as.numeric(names(vlist))
        pvs = sapply(vlist, function(v){
            pv = poly.pval2(y=y, poly=poly.nonfudged, v=v, sigma=sigma,bits=bits)$pv
        })

        return(data.frame(pvs=pvs,
                          locs=locs))
    }

    if(type=="cbs.noisy"){

        ## Draw new noise
        new.noise = rnorm(n,0,sigma.add)

        ## Get fudged sbs model
        stopifnot(numSteps%%2==0)
        h.fudged = circularBinSeg_fixedSteps(y + new.noise, numSteps=numSteps/2)
        poly.fudged = polyhedra(h.fudged)

        ## Get randomized p-value
        vlist <- make_all_segment_contrasts(h.fudged)
        vlist <- filter_vlist(vlist, visc)
        locs = as.numeric(names(vlist))
        pvs = sapply(vlist, function(v){
            pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                    sigma.add=sigma.add, orig.fudged.poly= poly.fudged,
                                    bits=bits, max.numIS=max.numIS, verbose=verbose,
                                    min.num.things=min.num.things)$pv})
        return(data.frame(pvs=pvs,
                          locs=locs))
    }

    if(type=="cbs"){
        ## Fit cbinseg on nonfudged data
        stopifnot(numSteps%%2==0)
        h.nonfudged = circularBinSeg_fixedSteps(y, numSteps=numSteps/2)
        poly.nonfudged = polyhedra(h.nonfudged)

        ## Get nonrandomized p-value
        vlist <- make_all_segment_contrasts(h.nonfudged)
        vlist <- filter_vlist(vlist, visc)
        locs = as.numeric(names(vlist))
        pvs = sapply(vlist, function(v){
        pv = poly.pval2(y=y, poly=poly.nonfudged, v=v, sigma=sigma, bits=bits)$pv
        })

        return(data.frame(pvs=pvs,
                          locs=locs))
    }
}

##' Does a single randomized wbs (rwbs) inference for a given y
##' @param numIntervals Defaults to \code{length(y)}, or if |intervals| is
##'     provided, the length of that.
##' @param intervals An object of class |intervals|.
##' @return List of simulation results (p-values, locations tested, etc.)
do_rwbs_inference <- function(y=y, max.numSteps=10, numIntervals=length(y),
                              intervals=NULL, consec=2,
                              sigma, postprocess=TRUE, how.close = 5,
                              better.segment=FALSE,
                              locs=1:length(y), numIS=100,
                              inference.type=inference.type,
                              improve.nomass.problem=TRUE, bits=1000,
                              max.numIS=2000,
                              write.time = FALSE, verbose=FALSE, mc.cores=1){

    ## Basic checks
    if(!is.null(intervals)){
        numIntervals = intervals$numIntervals
    }

    ## Fit initial WBS for a generous number of steps
    max.numSteps = 20
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

    ## ## Temporarily added
    ## if(length(retain)!=0) browser()

    ## Calculate the p-values
    vlist = vlist[retain]
    pvs = sapply(1:length(vlist), function(iv){
        printprogress(iv, length(vlist), type = "tests")
        v = vlist[[iv]]
        cumsum.v = cumsum(v)
        browser()
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
                                              mc.cores=mc.cores))
        if(write.time) write.time.to.file(myfile="rwbs-main-example-timing.txt")
        return(pv)})
    names(pvs) = names(vlist)
    return(list(pvs=pvs, locs.all=cp*cp.sign, locs.retained=as.numeric(names(pvs)), vlist=vlist) )
}


##' Does a single randomized wbs (rwbs) inference for a given y
do_rbs_inference <- function(y=y, max.numSteps=10, consec=2, sigma,
                             postprocess=TRUE, locs=1:length(y), numIS=100,
                             sigma.add = 0.2, bits=50, inference.type=c("rows", "pre-multiply"),
                             write.time=FALSE, numIntervals=length(y),
                             max.numIS=2000, verbose=FALSE){

    inference.type = match.arg(inference.type)

    ## Fit model and get IC information
    n = length(y)
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

    ## Retain only the changepoints we want results from:
    retain = which((cp %in% locs))
    if(length(retain)==0) return(list(pvs=c(), null.true=c()))
    vlist = vlist[retain]

    ## Do noise-added inference
    pvs = sapply(1:length(vlist), function(iv){
        v = vlist[[iv]]

        pv = randomize_addnoise(y= y, v=v, sigma=sigma, numIS=numIS,
                                sigma.add=sigma.add, orig.fudged.obj = h.fudged,
                                numSteps = stoptime+consec,
                                ic.poly = ic_obj$poly, bits=bits,
                                inference.type=inference.type,
                                max.numIS=max.numIS, verbose=verbose)

        if(write.time) write.time.to.file(myfile="rbs-main-example-timing.txt")
        return(pv)
    })
    names(pvs) = names(vlist)

    return(list(pvs=pvs, locs.all=cp*cp.sign, locs.retained=as.numeric(names(pvs))))
}


## do_rfl_inference <- function(y=y, max.numSteps=10, consec=2, sigma,
##                              postprocess=TRUE, locs=1:length(y), numIS=100,
##                              sigma.add = 0.2, bits=50){
##         ## Get nonrandomized p-value
##         D = genlassoinf::makeDmat(n,type='tf',ord=0)
##         f.nonfudged = genlassoinf::dualpathSvd2(y, D=D, maxsteps=1, approx=T)
##         Gobj.nonfudged = genlassoinf::getGammat.naive(obj=f.nonfudged, y=y, condition.step=1)
##         poly.nonfudged = polyhedra(obj=Gobj.nonfudged$G, u=Gobj.nonfudged$u)
##         vlist <- make_all_segment_contrasts(f.nonfudged)
##         if(!is.null(visc)){
##             retain = which((f.nonfudged$cp %in% visc))
##             if(length(retain)==0){
##                 return(data.frame(pvs=NA, locs=NA))
##             }
##             vlist = vlist[retain]
##         }
##         locs = (f.nonfudged$cp * f.nonfudged$cp.sign)[retain]

##         pvs = sapply(vlist, function(v){
##             pv = poly.pval2(y=y, poly=poly.nonfudged, v=v, sigma=sigma, bits=bits)$pv
##         })


##' Quick helper to write to file.
write.time.to.file <- function(myfile){
    line = Sys.time()
    write(toString(line),file=myfile,append=TRUE)
}

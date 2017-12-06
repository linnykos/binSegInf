## Synopsis: Compare FL and BS (nonrand) version to see why this doesn't work.



## Generate data
source("../main/wbs-tests/sim-helpers.R")
## type = match.arg(type)


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

nsim=400
results = mclapply(1:nsim,function(isim){printprogress(isim,nsim);dosim_compare_fl_and_bs()}, mc.cores=4)

## Detection of FL seems to be better than that of BS.
detection.sbs <- sapply(results, function(myresult){    return(length(myresult$pvs.sbs.nonrand)/4) })
detection.fl <- sapply(results, function(myresult){    return(length(myresult$pvs.fl.nonrand)/4) })
(mean(detection.sbs))
(mean(detection.fl))

## Conditional power
cond.pow.sbs <- sapply(results, function(myresult){return(mean(myresult$pvs.sbs.nonrand < 0.05/4)) })
cond.pow.fl <- sapply(results, function(myresult){return(mean(myresult$pvs.fl.nonrand < 0.05/4))})
cond.pow.fl <- cond.pow.fl[which(!is.nan(cond.pow.fl))]




## Fixing four steps might not be a great option for FL.
## Does this match what we saw with the generalized lasso paper?

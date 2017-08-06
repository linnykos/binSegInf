## Source in helpers
source("../main/justin/sim-helper.R")

## Takes a given mean and multiply the maximum to have noise*lev maximum
## height.
coriell_mn <- function(lev=1,n, std.bootstrap=1){
    ## newmn2 = (newmn[1101:1300][seq(from=1,to=200,length=100)])
    h = max(abs(newmn))
    return((newmn / h * std.bootstrap) * lev)
}


onesim_bsft <- function(sim.settings){

    ## Reassign things
    sigma = sim.settings$sigma
    lev = sim.settings$lev
    nsim.is = sim.settings$nsim.is
    numSteps = sim.settings$numSteps
    numIntervals = sim.settings$numIntervals
    n = sim.settings$n
    meanfun = sim.settings$meanfun
    reduce = sim.settings$reduce
    augment = sim.settings$augment
    bootstrap = sim.settings$bootstrap
    std.bootstrap = sim.settings$std.bootstrap
    cleanmn.bootstrap = sim.settings$cleanmn.bootstrap
    thresh = sim.settings$thresh

    ## Generate data (same across all onesim_OOO functions)
    my.mn <- meanfun(lev,n)
    y = my.mn + stats::rnorm(n,0,sigma)

    ## Do SBS-FS inference (basic structure is similar)
    method <- binSeg_fixedThresh
    obj <- method(y, thresh=thresh)
    poly <- polyhedra(obj)
    if(length(obj$cp)==0){
        return(NULL)
    } else{
        contrasts <- make_all_segment_contrasts(obj)
        pvec = sapply(1:length(obj$cp), function(ii){
            poly.pval(y= y, G= poly$gamma, u=poly$u, v=contrasts[[ii]],sigma=sigma,
                      bits=100)$pv
        })
        names(pvec) = obj$cp
        return(pvec)
    }
}

onesim_bsfs <- function(sim.settings){

    ## Reassign things
    sigma = sim.settings$sigma
    lev = sim.settings$lev
    nsim.is = sim.settings$nsim.is
    numSteps = sim.settings$numSteps
    numIntervals = sim.settings$numIntervals
    n = sim.settings$n
    meanfun = sim.settings$meanfun
    reduce = sim.settings$reduce
    augment = sim.settings$augment

    ## Generate data (same across all onesim_OOO functions)
    my.mn <- meanfun(lev,n)
    y = (my.mn + stats::rnorm(n,0,sigma))

    ## Do SBS-FS inference (basic structure is similar)
    method <- binSeg_fixedSteps
    obj <- method(y, numSteps)
    poly <- polyhedra(obj)
    contrasts <- make_all_segment_contrasts(obj)
    pvec = sapply(1:length(obj$cp), function(ii){
        poly.pval(y= y, G= poly$gamma, u=poly$u, v=contrasts[[ii]],sigma=sigma,
                  bits=100)$pv
    })
    names(pvec) = obj$cp
    return(pvec)
}


onesim_wbs <- function(sim.settings){

    ## Reassign things
    sigma = sim.settings$sigma
    lev = sim.settings$lev
    nsim.is = sim.settings$nsim.is
    numSteps = sim.settings$numSteps
    numIntervals = sim.settings$numIntervals
    n = sim.settings$n
    meanfun = sim.settings$meanfun
    reduce = sim.settings$reduce
    augment = sim.settings$augment
    bootstrap = sim.settings$bootstrap
    std.bootstrap = sim.settings$std.bootstrap
    cleanmn.bootstrap = sim.settings$cleanmn.bootstrap
    type = sim.settings$type

    ## sigma=1; lev=0; nsim.is=100; numSteps=1;
    ## numIntervals=20; n=6; meanfun=onejump;
    ## reduce=TRUE;augment=TRUE;  bootstrap=FALSE; std.bootstrap=NULL;
    ## cleanmn.bootstrap=NULL;
    ## type = "randomized"


    ## Generate data
    y = meanfun(lev,n) + stats::rnorm(n,0,sigma)

    ## Do WBS_FS inference (basic structure is similar)
    method <- wildBinSeg_fixedSteps
    intervals <- generate_intervals(length(y),numIntervals)
    obj <- method(y, numSteps=numSteps, intervals=intervals)
    contrasts <- make_all_segment_contrasts(obj)
    pvec = pvec.plain = setNames(rep(NA,length(obj$cp)), obj$cp)
    for(ii in 1:length(obj$cp)){
        poly <- polyhedra(obj, v = contrasts[[ii]], reduce=reduce, sigma=sigma)
        if(type=="plain"){
            pvec.plain[ii] <- poly.pval2(y=y, poly=poly, v=contrasts[[ii]],
                                         sigma=sigma, reduce=reduce)$pv
        } else {
            pvec[ii] <- randomized_wildBinSeg_pv(y=y,
                                                 v=contrasts[[ii]], sigma=sigma,
                                                 numSteps=numSteps,
                                                 numIntervals=numIntervals,
                                                 nsim.is=nsim.is, bits=100,
                                                 reduce=reduce,
                                                 augment=augment)
        }
    }
    if(type=="plain"){ return(pvec.plain) } else {  return(pvec) }
}



onesim_naive <- function(sim.settings){

    ## Reassign things
    sigma = sim.settings$sigma
    lev = sim.settings$lev
    nsim.is = sim.settings$nsim.is
    numSteps = sim.settings$numSteps
    numIntervals = sim.settings$numIntervals
    n = sim.settings$n
    meanfun = sim.settings$meanfun
    reduce = sim.settings$reduce
    augment = sim.settings$augment
    bootstrap = sim.settings$bootstrap
    std.bootstrap = sim.settings$std.bootstrap
    cleanmn.bootstrap = sim.settings$cleanmn.bootstrap
    type = sim.settings$type

    ## Generate data
    y = meanfun(lev,n) + stats::rnorm(n,0,sigma)

    ## Do binseg naive z-tests
    method <- binSeg_fixedSteps
    obj <- method(y, numSteps)
    pvec.naive = sapply(1:length(obj$cp), function(ii){
        ztest(contrast_vector(obj, ii), y, sigma, 0.05)})
    return(pvec.naive)
}



onesim_fusedlasso <- function(sim.settings){

    ## Reassign things
    sigma = sim.settings$sigma
    lev = sim.settings$lev
    nsim.is = sim.settings$nsim.is
    numSteps = sim.settings$numSteps
    numIntervals = sim.settings$numIntervals
    n = sim.settings$n
    meanfun = sim.settings$meanfun
    reduce = sim.settings$reduce
    augment = sim.settings$augment
    bootstrap = sim.settings$bootstrap
    std.bootstrap = sim.settings$std.bootstrap
    cleanmn.bootstrap = sim.settings$cleanmn.bootstrap
    type = sim.settings$type

    ## Generate data (same across all onesim_OOO functions)
    y = meanfun(lev,n) + stats::rnorm(n,0,sigma)

    ## Do fused lasso inference (basic structure is similar)
    D = genlassoinf::makeDmat(n,type='tf',ord=0)
    obj <- genlassoinf::dualpathSvd2(y,D=D,maxsteps=numSteps,approx=TRUE)
    contrasts <- make_all_segment_contrasts(obj)
    pvec = pvec.plain = setNames(rep(NA,length(obj$cp)), obj$cp)
    poly <- polyhedra(obj$Gobj.naive$G, obj$Gobj.naive$u)
    for(ii in 1:length(obj$cp)){
       if(type=="plain"){
            pvec.plain[ii] <- poly.pval2(y=y, poly=poly, v=contrasts[[ii]],
                                         sigma=sigma, reduce=reduce)$pv
        } else {
            pvec[ii] <- randomized_genlasso_pv(y=y, v=contrasts[[ii]],
                                               sigma=sigma, numSteps=numSteps,
                                               numIntervals=numIntervals,
                                               nsim.is=nsim.is, bits=100,
                                               reduce=reduce, augment=augment)
        }
    }
    if(type=="plain"){
        return(pvec.plain)
    } else {
        return(pvec)
    }
}



## ##' Simulation driver.
## ##' @param sim.setting list of simulation settings, set externally.
## ##' @param filename name of R data file to save in.
## ##' @import parallel
## sim_driver <- function(sim.settings, filename, dir="../data",seed=NULL,
##                        mc.cores=4, reduce=FALSE, resid.cleanmn=NULL){
##     levs = sim.settings$levs
##     n.levs = length(levs)
##     results <- replicate(n.levs, list())
##     names(results) = paste("jump size",levs)
##     for(i.lev in 1:n.levs){

##         ## Run simulations
##         cat("signal strength (level)", i.lev, "out of", n.levs, fill=TRUE)
##         nsim = sim.settings$nsims[i.lev]

##         ## Add tryCatch()
##         manysimresult =
##            ## mclapply(1:nsim, function(isim){
##            lapply(1:nsim, function(isim){
##                cat("\r", "simulation ", isim, "out of", nsim)
##                ## cat("simulation ", isim, "out of", nsim)
##                onesim(isim, lev=sim.settings$levs[i.lev],
##                       sigma=sim.settings$sigma,
##                       nsim.is=sim.settings$nsim.is,
##                       numSteps=sim.settings$numSteps,
##                       numIntervals=sim.settings$numIntervals,
##                       n=sim.settings$n,
##                       mn=sim.settings$mn,
##                       seed=isim,
##                       reduce=reduce,
##                       bootstrap=sim.settings$bootstrap,
##                       std.bootstrap=sim.settings$std.bootstrap,
##                       augment = sim.settings$augment,
##                       resid.cleanmn = resid.cleanmn
##                       )
##                })
##            ## }, mc.cores = mc.cores)

##         ## Extract plist
##         plist.bsfs <- lapply(manysimresult, function(a)a$p.bsfs)
##         plist.wbsfs <- lapply(manysimresult, function(a)a$p.wbsfs)
##         plist.wbsfs.plain <- lapply(manysimresult, function(a)a$p.wbsfs.plain)

##         ## Reformat to pmat
##         pmat.bsfs = reformat(plist.bsfs, n=sim.settings$n)
##         pmat.wbsfs = reformat(plist.wbsfs, n=sim.settings$n)
##         pmat.wbsfs.plain = reformat(plist.wbsfs.plain, n=sim.settings$n)

##         ## Store results
##         results[[i.lev]] <- list(pmat.bsfs = pmat.bsfs,
##                                  pmat.wbsfs =  pmat.wbsfs,
##                                  pmat.wbsfs.plain = pmat.wbsfs.plain)

##         save(results, sim.settings, file = file.path(dir,filename))
##         cat(fill=TRUE)
##     }
## }

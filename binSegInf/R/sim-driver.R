##' Generates one-jump mean
mn.onejump <- function(lev,n){c(rep(0,n/2),rep(lev,n/2))}
##' Generates two-jump mean
mn.twojump <- function(lev,n){c(rep(0,n/3),rep(lev,n/3), rep(0,n/3))}

##' Simulation inner function.
onesim <- function(isim, sigma, lev, nsim.is, numSteps, numIntervals, n, mn,
                   seed=NULL,reduce, bootstrap=FALSE, std=NULL, augment, resid.cleanmn){

    ## generate data
    if(!is.null(seed)) set.seed(seed)
    my.mn <- mn(lev,n)
    if(is.null(bootstrap)) bootstrap = FALSE
    if(bootstrap) y = (my.mn + bootstrap_sample(resid.cleanmn, seed=seed))
    if(!bootstrap) y = (my.mn + stats::rnorm(n,0,sigma))

    ###########################
    ## Do SBS-FS inference ####
    ###########################
    method <- binSeg_fixedSteps
    obj <- method(y, numSteps)
    poly <- polyhedra(obj)
    p.bsfs = rep(NA,length(obj$cp))
    names(p.bsfs) = obj$cp
    contrast = list()
    for(ii in 1:length(obj$cp)){
        contrast[[ii]] <- contrast_vector(obj, ii)
        p.bsfs[ii] <- poly.pval(y= y, G= poly$ gamma, u=poly$u,
                                v=contrast[[ii]],sigma=sigma, bits=100)$pv
    }

    ###########################
    ## Do WBS-FS inference ####
    ###########################
    method <- wildBinSeg_fixedSteps
    intervals <- generate_intervals(n=length(y),numIntervals=numIntervals,seed=seed)
    obj <- method(y, numSteps=numSteps, intervals=intervals)
    contrast <- make_all_segment_contrasts(obj)
    p.wbsfs = p.wbsfs.plain = rep(NA,length(obj$cp))
    for(ii in 1:length(obj$cp)){
        system.time(poly <- polyhedra(obj, v = contrast[[ii]], reduce=TRUE, sigma=sigma))
        p.wbsfs.plain[ii] <- poly.pval2(y=y, poly=poly, v=contrast[[ii]],
                                        sigma=sigma, reduce=reduce)$pv
        p.wbsfs[ii] <- randomized_wildBinSeg_pv(y=y,
                                                v=contrast[[ii]], sigma=sigma,
                                                numSteps=numSteps,
                                                numIntervals=numIntervals,
                                                nsim.is=nsim.is, bits=100,
                                                reduce=reduce,
                                                augment=augment)
    }
    names(p.bsfs) = names(p.wbsfs.plain) = names(p.wbsfs) = obj$cp

    ########################
    ## Do CBS inference ####
    ########################

    return(list(p.wbsfs.plain = p.wbsfs.plain,
                p.bsfs=p.bsfs,
                p.wbsfs=p.wbsfs))
}


##' Simulation driver.
##' @param sim.setting list of simulation settings, set externally.
##' @param filename name of R data file to save in.
##' @import parallel
sim_driver <- function(sim.settings, filename, dir="../data",seed=NULL,
                       mc.cores=4, reduce=FALSE, resid.cleanmn=NULL){
    levs = sim.settings$levs
    n.levs = length(levs)
    results <- replicate(n.levs, list())
    names(results) = paste("jump size",levs)
    for(i.lev in 1:n.levs){

        ## Run simulations
        cat("signal strength (level)", i.lev, "out of", n.levs, fill=TRUE)
        nsim = sim.settings$nsims[i.lev]

        manysimresult =
           mclapply(1:nsim, function(isim){
               cat("\r", "simulation ", isim, "out of", nsim)
               ## cat("simulation ", isim, "out of", nsim)
               onesim(isim, lev=sim.settings$levs[i.lev],
                      sigma=sim.settings$sigma,
                      nsim.is=sim.settings$nsim.is,
                      numSteps=sim.settings$numSteps,
                      numIntervals=sim.settings$numIntervals,
                      n=sim.settings$n,
                      mn=sim.settings$mn,
                      seed=isim,
                      reduce=reduce,
                      bootstrap=sim.settings$bootstrap,
                      std=sim.settings$std,
                      augment = sim.settings$augment,
                      resid.cleanmn = resid.cleanmn
                      )
           }, mc.cores = mc.cores)

        ## Extract plist
        plist.bsfs <- lapply(manysimresult, function(a)a$p.bsfs)
        plist.wbsfs <- lapply(manysimresult, function(a)a$p.wbsfs)
        plist.wbsfs.plain <- lapply(manysimresult, function(a)a$p.wbsfs.plain)

        ## Reformat to pmat
        pmat.bsfs = reformat(plist.bsfs, n=sim.settings$n)
        pmat.wbsfs = reformat(plist.wbsfs, n=sim.settings$n)
        pmat.wbsfs.plain = reformat(plist.wbsfs.plain, n=sim.settings$n)

        ## Store results
        results[[i.lev]] <- list(pmat.bsfs = pmat.bsfs,
                                 pmat.wbsfs =  pmat.wbsfs,
                                 pmat.wbsfs.plain = pmat.wbsfs.plain)

        save(results, sim.settings, file = file.path(dir,filename))
        cat(fill=TRUE)
    }
}

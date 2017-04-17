##' Toy data generating function, of up-then-down two-jump signal, for use in
##' Rmd simulation files.
##' @param seed Seed to generate data from
##' @param n Data length
##' @param lev jump size
##' @param sigma Data noise
##' @param type data
twojump_data <- function(seed=NULL,n,lev,sigma, type=c("data","mean")){

    type = match.arg(type)

    ## Basic checks
    if(n%%3!=0) stop("Provide |n| divisible by 3 !")

    ## Generate mean
    mn <- rep(c(0,lev,0), each=n/3)
    if(type=="mean")  return(mn)

    ## Generate data
    if(!is.null(seed))  set.seed(seed)
    y <- mn + rnorm(n, 0, sigma)
    if(type=="data")  return(y)

}


##' Toy data generating function, of upward one-jump signal, for use in Rmd
##'  simulation files.
##' @param seed Seed to generate data from
##' @param n Data length
##' @param lev jump size
##' @param sigma Data noise
##' @param type data
onejump_data <- function(seed=NULL, n, lev,sigma, type=c("data","mean")){

   type = match.arg(type)

    ## Basic checks
    if(n%%2!=0) stop("Provide |n| divisible by 3 !")

    ## Generate mean
    mn <- rep(c(0,lev), each=n/2)
    if(type=="mean")  return(mn)

    ## Generate data
    if(!is.null(seed))  set.seed(seed)
    y <- mn + rnorm(n, 0, sigma)
    if(type=="data")  return(y)

}







##' See if the vector of p-values (whose rownames are locations) \code{pv} has
##' location \code{loc}.
##' @param pv numeric vector of p-values
##' @param loc location
has.loc <- function(pv,loc){loc %in% names(pv)}

##' Create verdicts from the vector of p-values (whose rownames are locations)
##' \code{pv} has location \code{loc}, accounting for how many tests were
##' conducted.
##' @param pv numeric vector of p-values
##' @param loc location
rejected <- function(pv){pv<(0.05/length(pv))}


##' Wrapper to pvalue()
pvalue2 = function(v,...){pvalue(contrast=v,...)}


##' Get power from pvalues
##' @param pvals list of numeric vector of p-values
##' @param loc location to condition on (can be a vector)
get.power <- function(pvals, loc=n/2){

    ## preprocess
    pvals = lapply(pvals, function(mypvals){if(class(mypvals)=="list") unlist(mypvals) else mypvals})

    ## do stuff
    verdicts = do.call(c,lapply(pvals, rejected))
    ind = (names(verdicts) %in% loc)
    pow = sum(verdicts[ind])/length(ind)
    return(pow)
}


##' Helper to inquire about changepoint being tested, \code{cp}, being contained
##' in the interval set $i$.
##' @param i changepoint set
##' @param cp single location
##' @return TRUE if \code{cp} is covered by \code{i}
i.covers.cp <- function(i,cp){
    contained = (i$starts <= cp & cp<i$ends)
    return(any(contained))
}



##' filter NULLs out of a pvmat list.
.filternull <- function(pvmat){
    emptyguys = lapply(pvmat, function(pvobj) return(length(pvobj)==0))
    return(pvmat[emptyguys])
}


## ##' Single simulation driver, for experiments. Works for wildBinSeg, binSeg_OOO.
## one_sim <- function(pfun, method, lev, datfun, seed=NULL,n, sigma,
##                     thresh=NULL,size=NULL,pfun2=NULL, ...){
    
##     ## Generate data
##     y = datfun(seed,n,lev,sigma)
    
##     ## Run method /once/, collect things
##     obj = method(y,thresh, numIntervals=numIntervals, verbose=FALSE)
##     if(length(obj$cp)==0 | all(is.na(obj$cp)))return(NULL)
    
##     ## Collect a polyhedron
##     poly <- polyhedra(obj)
    
##     ## Make contrasts
##     vs = make_all_segment_contrasts(obj)

##     ## Do inference and return
##     pv = sapply(vs, function(v){
##         pfun(y=y,v=v,poly=poly,sigma=sigma,...)})

##     ## Optionally, compare it directly to a different pfun2
##     if(!is.null(pfun2)){
##         pv2 = sapply(vs, function(v){ pfun2(y=y,v=v,poly=poly,sigma=sigma,...)})
##         return(c(pv,pv2))
##     }
    
##     return(pv)
## }


##' Plotter for list of pval vectors.
plot_pvals <- function(pvals.list, title=c("wbs","sbs")){
    invisible(lapply(pvals.list, function(plist){
        ii <<- ii+1
        
        ## Extract verdicts
        pp = unlist(plist)
        vd = unlist(lapply(plist, function(pvec){pvec<0.05/length(pvec)}))
        
        ## Get indices of approximately recovered locations
        prox = 0
        ind = which(names(pp) %in% c(n/3 + (-prox):prox, 2*n/3 + (-prox):prox))
        
        ## Make qqplot
        qqunif(unlist(pp[ind]))
        title(paste("Jump size =", levs[ii]))

        ## Calculate power
        title(sub=round(sum(vd[ind])/length(ind),3))
        
    }))}

##' Generates one-jump mean
mn.onejump <- function(lev,n){c(rep(0,n/2),rep(lev,n/2))}
##' Generates two-jump mean
mn.twojump <- function(lev,n){c(rep(0,n/3),rep(lev,n/3), rep(0,n/3))}

##' Simulation inner function.
onesim <- function(isim, sigma, lev, nsim.is, numSteps, numIntervals, n, mn){

    ## generate data
    y <- mn(lev,n) + rnorm(n,0,sigma)
    set.seed(0)
    y <-runif(10)

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
        p.bsfs[ii] <- poly.pval(y=y,
                                G=poly$ gamma,u=poly$u, v=contrast[[ii]],sigma=sigma, bits=100)$pv
    }
    p.bsfs = cbind(rep(isim,length(obj$cp)), obj$cp, p.bsfs)

    ###########################
    ## Do WBS-FS inference ####
    ###########################
    method <- wildBinSeg_fixedSteps
    intervals <- generate_intervals(n=length(y),numIntervals=numIntervals)
    obj <- method(y, numSteps=numSteps, intervals=intervals)
    poly <- polyhedra(obj)
    p.wbsfs = p.wbsfs.nonrand = rep(NA,length(obj$cp))
    contrast <- make_all_segment_contrasts(obj)
    for(ii in 1:length(obj$cp)){
        p.wbsfs.nonrand[ii] <- poly.pval(y=y, G=poly$gamma, u=poly$u,
                                         v=contrast[[ii]], sigma=sigma)$pv
        p.wbsfs[ii] <- randomized_wildBinSeg_pv(y=y,
                                                v=contrast[[ii]], sigma=sigma,
                                                numSteps=numSteps,
                                                numIntervals=numIntervals,
                                                nsim.is=nsim.is, bits=100)
    }
    p.wbsfs = cbind(rep(isim,length(obj$cp)), obj$cp, p.wbsfs)
    p.wbsfs.nonrand = cbind(rep(isim,length(obj$cp)), obj$cp, p.wbsfs.nonrand)
    colnames(p.bsfs) = colnames(p.wbsfs.nonrand) = colnames(p.wbsfs) = c("isim","cp","pv")
    return(list(p.wbsfs.nonrand = p.wbsfs.nonrand,
                p.bsfs=p.bsfs,
                p.wbsfs=p.wbsfs))
}


##' Simulation driver.
sim_driver <- function(sim.settings, filename, dir="../data"){
    levs = sim.settings$levs
    n.levs = length(levs)
    results <- replicate(n.levs, list())
    names(results) = paste("jump size",levs)
    for(i.lev in 1:n.levs){
        ## Run simulations
        cat(i.lev, "out of", n.levs, fill=TRUE)
        manysimresult = lapply(1:sim.settings$nsims[i.lev],
                               function(isim){
                                   if(isim%%100==1){print(isim)}
                                   onesim(isim, lev=sim.settings$levs[i.lev],
                                          sigma=sim.settings$sigma,
                                          nsim.is=sim.settings$nsim.is,
                                          numSteps=sim.settings$numSteps,
                                          numIntervals=sim.settings$numIntervals,
                                          n=sim.settings$n,
                                          mn=sim.settings$mn)})
        ## Extract and Aggregate
        browser()
        plist.bsfs <- lapply(manysimresult, function(a)a$p.bsfs)
        plist.wbsfs <- lapply(manysimresult, function(a)a$p.wbsfs)
        plist.wbsfs.nonrand <- lapply(manysimresult, function(a)a$p.wbsfs.nonrand)
        pmat.bsfs = do.call(rbind, plist.bsfs)
        pmat.wbsfs = do.call(rbind, plist.wbsfs)
        pmat.wbsfs.nonrand = do.call(rbind, plist.wbsfs.nonrand)
        results[[i.lev]] <- list(pmat.bsfs = pmat.bsfs,
                                 pmat.wbsfs =  pmat.wbsfs,
                                 pmat.wbsfs.nonrand = pmat.wbsfs.nonrand)

        save(results, sim.settings, file = file.path(dir,filename))
    }
}

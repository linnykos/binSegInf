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
onesim <- function(isim, sigma, lev, nsim.is, numSteps, numIntervals, n, mn,
                   seed=NULL,reduce, bootstrap=FALSE, std=NULL){

    ## generate data
    if(!is.null(seed)) set.seed(seed)
    my.mn <- mn(lev,n)
    y <- my.mn + rnorm(n,0,sigma)

    if(is.null(bootstrap)) bootstrap = FALSE
    if(bootstrap) y = (my.mn + bootstrap_sample(resid.cleanmn, seed=seed))
    if(!bootstrap) y = (my.mn + rnorm(length(y),0,sigma))

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
    ## p.bsfs = cbind(rep(isim,length(obj$cp)), obj$cp, p.bsfs)

    ###########################
    ## Do WBS-FS inference ####
    ###########################
    method <- wildBinSeg_fixedSteps
    intervals <- generate_intervals(n=length(y),numIntervals=numIntervals)
    obj <- method(y, numSteps=numSteps, intervals=intervals)
    contrast <- make_all_segment_contrasts(obj)
    p.wbsfs = p.wbsfs.plain = rep(NA,length(obj$cp))
    for(ii in 1:length(obj$cp)){
        poly <- polyhedra(obj, v = contrast[[ii]], reduce=reduce)
        p.wbsfs.plain[ii] <- poly.pval(y=y, G=poly$gamma, u=poly$u,
                                         v=contrast[[ii]], sigma=sigma)$pv
        p.wbsfs[ii] <- randomized_wildBinSeg_pv(y=y,
                                                v=contrast[[ii]], sigma=sigma,
                                                numSteps=numSteps,
                                                numIntervals=numIntervals,
                                                nsim.is=nsim.is, bits=100,
                                                reduce=reduce)
    }
    ## p.wbsfs = cbind(rep(isim,length(obj$cp)), obj$cp, p.wbsfs)
    ## p.wbsfs.plain = cbind(rep(isim,length(obj$cp)), obj$cp, p.wbsfs.plain)
    ## colnames(p.bsfs) = colnames(p.wbsfs.plain) = colnames(p.wbsfs) = c("isim","cp","pv")
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
                       mc.cores=4, reduce=FALSE){
    levs = sim.settings$levs
    n.levs = length(levs)
    results <- replicate(n.levs, list())
    names(results) = paste("jump size",levs)
    ## ptm <- proc.time()
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
                      seed=seed,
                      reduce=reduce,
                      bootstrap=sim.settings$bootstrap,
                      std=sim.settings$std
                      )
               ## print(proc.time() - ptm)
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


##' Gets sigma by fitting basic wbs::sbs on it, with default settings for
##' complexity.
##' @param y data vector
##' @import wbs
##' @export
get_sigma <- function(y){
  cps = sort(changepoints(wbs::sbs(y))$cpt.th[[1]])
  segment.inds = sapply(1:(length(cps)+1),
                      function(ii) (c(0,cps,length(y))[ii]+1):(c(0,cps,length(y))[ii+1]))
  mn = rep(NA,length(y))
  for(ind in segment.inds) mn[ind] <- mean(y[ind])
  return(sd(y-mn))
}

##' Return piecewise mean.
##' @param y data vector.
##' @param cp changepoint vector. Assumed to be sorted.
##' @import glmgen
##' @export
get_piecewise_mean <- function(y, cp){
  if(!all.equal(sort(cp), cp)) stop ("cp is not sorted!")
  aug.cp = c(0,cp,length(y))
  segment.inds = sapply(1:(length(cp)+1),
                      function(ii){ (aug.cp[ii]+1):aug.cp[ii+1]})
  mn = rep(NA,length(y))
  for(ind in segment.inds) mn[ind] <- mean(y[ind])
  return(mn)
}



##' Naive inference
##' @param v contrast
##' @param y data
##' @param sigma noise level
##' @return p-value for one-sided test of \eqn{H_0: v^Ty=0}.
ztest <- function(v,y,sigma,alpha=0.05){
  ## stopifnot(sum(v*y)>0)
  if(sum(v*y)<=0){print(sum(v*y));print("v'y is not positive!!!! Something is slightly fishy.")}
  test.stat <- sum(v*y) * 1/sqrt(sum(v*v)) * 1/sigma
  ## @param alpha significance level
  ## cutoff <- qnorm(1-alpha, 0,1)
  return(dnorm(test.stat))
}



##' Reformat then plot the p-values, from plist
reformat <- function(my.plist, n){
    my.pmat = do.call(Matrix::rBind, lapply(my.plist,
                                    function(x){
                                      myrow <- rep(NA,n);
                                      myrow[as.numeric(names(x))] <- x
                                      myrow
                                      }))

## lapply(plist.bsfs, function(x){
##     myrow <- rep(NA,n);
##     myrow[as.numeric(names(x))] <- x
##     myrow
## })
  colnames(my.pmat) = 1:n
  return(my.pmat)
}

##' Get the p-values at locations, from pmat
get_condit_pvals <- function(my.pmat, loc){
  pvals = as.numeric(my.pmat[,loc])
  return(pvals)
}

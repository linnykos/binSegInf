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
get.power <- function(pvals, loc){

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
  return(stats::sd(y-mn))
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
  return(stats::dnorm(test.stat))
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


## Generates one/two-jumped means
onejump <- function(lev,n){c(rep(0,n/2),rep(lev,n/2))}
twojump <- function(lev,n){c(rep(0,n/3),rep(lev,n/3), rep(0,n/3))}

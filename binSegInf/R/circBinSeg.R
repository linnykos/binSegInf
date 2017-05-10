##' Algorithm for Circular binarysegmentation
##' @param y Data vector.
circBinSeg <- function(y, return.env=FALSE){

    ## Initialize things
    n = length(y)
    env = new.env()
    env$cplist2 = cplist2(n)
    env$y = y

    ## Run CBS
    circBinSeg_inner(y=y, s=1, e=n, j=0,k=1, n=n,env=env)

    ## Gather output from |env| and return it.
    bs.output = list(cplist2 = cplist2,
                     y = y)

    obj = structure(list(bs.output = bs.output,
                         y = y),
                         ## cp = trim.vec(env$blist$mat[,3]),
                         ## cp.sign = trim.vec(env$zlist$mat[,3])) ,
                    class = "cbs")

    if(return.env){ return(env) } else{ return(obj) }
}



##' Inner function for binary segmentation with fixed threshold. The wrapper
##' \code{circBinSeg()} is intended to be used by user.
##' @param s Starting index, in the vector-valued data. Must be an integer
##'     larger than or equal to 1, and strictly smaller than \code{e}.
##' @param e Ending index, in the vector-valued data. Must be an integer smaller
##'     than or equal to \code{n}, and strictly larger than \code{s}.
##' @param j The depth of the recursion on hand.
##' @param k The indexing of the node location, from left to right, in the
##'     \emph{complete} binary tree.
##' @param y The original data.
##' @param n The length of the data \code{y}.
##' @param env environment where the algorithm output lives.
##' @param verbose \code{TRUE} if you want this to be loud.
circBinSeg_inner <- function(y, s, e, j, k, n, env,verbose=FALSE){

    ## If segment is 2-lengthed, terminate
    if(e-s<1){
        return(env )

    ## Otherwise, calculate CUSUMs
    } else {
        ## Obtain breakpoints and direction
        cobj <- cbs_maximize(env$y, s, e)
        snew <- cobj$snew
        enew <- cobj$enew
        z <- sign(cobj$maxcrit)
        ## sigma=1
        ## thresh <- makethresh(1,s,e,n,sigma)
        thresh=0

        ## Check threshold exceedance, then store
        if(abs(cobj$maxcrit) < thresh){
            addrow(env$cplist2, j=j, k=k, s=snew, e=enew, s0=s, e0=e, z=z, pass=FALSE)
            return(env)
        } else {

        if(verbose) cat("Changepoints", snew, enew,
                        "passed!, between s and e:",s, e, fill=T)

            addrow(env$cplist2, j=j, k=k, s=snew, e=enew, s0=s, e0=e, z=z, pass=TRUE)

            ## Recurse
            circBinSeg_inner(y=y,s=s,    e=snew, j=j+1, k=3*k,   env=env)
            circBinSeg_inner(y=y,s=snew+1, e=enew, j=j+1, k=3*k+1, env=env)
            circBinSeg_inner(y=y,s=enew+1, e=e,  j=j+1, k=3*k+2, env=env)
        }
    }
}


##' Computes the CBS CUSUM (hump) statistic.  Note, we calculate this as an
##' upward hump by default.
##' @param s starting index.
##' @param snew first breakpoint index.
##' @param enew second jbreakpoint index.
##' @param e end index.
##' @param y data.
##' @param n length of data. Defaults to \code{length(y)}.
##' @param upward Whether you want upward jump in the cusum
##'     calculation. Defaults to TRUE.
##' @param contrast.vec If TRUE, then the contrast vector v for cusum=v'y is
##'     returned.
##' @param unsigned if TRUE, then returns the contrast vector that makes the
##'     contrast with |y| positive, or the absolute value of the
##'     contrast. Defaults to FALSE.
##' @export
crit <- function(y,s=1,snew,enew,e=length(y), n=length(y), upward=TRUE, unsigned = FALSE, contrast.vec=FALSE){


    ## Basic checks
    if(enew>=e) stop("enew must be strictly smaller than e!")
    if(snew>=enew) stop("enew must be strictly larger than snew!")
    if(s>snew) stop("snew must be larger than or equal to s!")
    if(unsigned & is.null(y)) stop("Can't produce an unsigned linear contrast or contrast vector without a y")
    if(is.null(y)){ if(is.null(n)) stop("Must provide one of y or n.")}
    if(is.null(n)) n = length(y)

    ## Calculate things
    ind.middle = (snew+1):enew
    ind.outer = c(s:snew, (enew+1):e)
    len.middle = length(ind.middle)
    len.outer =  length(ind.outer)

    ## Make contrast vector
    v = rep(0,length(y))
    v[ind.middle] = 1/len.middle
    v[ind.outer]  = -1/len.outer
    v = v * sqrt(1/(1/len.middle + 1/len.outer))
    if(!upward) v = -v
    if(unsigned) v = v * sign(sum(v*y))

    ## Return the right thing
    if(contrast.vec){
        return(v)
    } else if(!is.null(y)) {
        return(sum(v*y))
    } else {
        stop("Check your options for cusum() again!")
    }
}


##' Do CBS maximization.
##' @param y data vector
##' @param s start index
##' @param e end index
cbs_maximize <- function(y, s, e){

    ## Calculate cusums
    ssmat <- t(utils::combn(s:(e-1),2))
  crits <- apply(ssmat, 1, function(ss){
    crit(y=y, s=s, snew=ss[1], enew=ss[2], e=e, unsigned = FALSE)})
    ind.max <- which.max(abs(crits))
    ss = ssmat[ind.max,]
    return(list(snew = ss[1], snew = ss[2], maxcrit = crits[ind.max]))
}


makethresh <- function(s,s0,e0,n,sigma){
  ## Make it so you only have to query from the pre-made table
  return(0)
}

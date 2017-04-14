##' Algorithm for Circular binary segmentation
##' @param y Data vector.
circBinSeg <- function(y){

    ## Initialize things 
    env = new.env()
    env$Blist = env$Zlist = env$passmat = cplist()
    env$y = y

    ## Run CBS
    circBinSeg_inner(y=y, s=1, e=n, env=env)

    ## Gather output from |env| and return it.
    bs.output = list(slist = trim(env$slist),
                     elist = trim(env$elist),
                    selist = trim(env$selist),
                     blist = trim(env$blist),
                     Blist = trim(env$Blist),
                     Zlist = trim(env$Zlist),
                     cp = trim.vec(env$blist$mat[,3]),
                     cp.sign = trim.vec(env$zlist$mat[,3]),
                     y = y)
    
    obj = structure(list(bs.output = bs.output,
                         y = y,
                         cp = trim.vec(env$blist$mat[,3]),
                         cp.sign = trim.vec(env$zlist$mat[,3])) ,
                    class = "cbs")
    
    if(return.env){ return(env) } else{ return(obj) }

}

##' Inner function for binary segmentation with fixed threshold. The wrapper
##' \code{circBinSeg()} is intended to be used by user.
circBinSeg_inner <- function(y, s, e, env){

    ## If segment is 2-lengthed, terminate
    if(e-s<1){
        return() 
    ## Otherwise, calculate CUSUMs
    } else {

        ## Obtain breakpoints and direction
        cobj <- cbs_maximize(y,s,e)
        s1 <- cobj$s1
        s2 <- cobj$s2
        z <- sign(cobj$maxcusumm) 
        if(verbose) cat("the (potential) changepoints", s1, s2,
                        "are being considered, between s and e:",s, e, fill=T)
        
        thresh <- makethresh(1,s,e,n,sigma)
        ## Check threshold exceedance, then store
        if(abs(all.cusums[b]) < thresh){

            env$Blist = add(env$Blist,j+1,k,b)
            env$Zlist = add(env$Zlist,j+1,k,z)
            env$slist = add(env$slist,j,  k,s)
            env$elist = add(env$elist,j,  k,e)
            env$selist = add(env$selist, j  ,k, paste(s,e))
            
            return(env)
        } else { 
            if(verbose) cat("the biggest cusum was", all.cusums[b],
                            "which passed the threshold:", thresh,
                            "frm between s and e:",s,e, fill=T)

            env$blist = add(env$blist, j+1,k, b)
            env$Blist = add(env$Blist, j+1,k, b)
            env$zlist = add(env$zlist, j+1,k, z)
            env$Zlist = add(env$Zlist, j+1,k, z)
            env$slist = add(env$slist, j  ,k, s)
            env$elist = add(env$elist, j  ,k, e)
            env$selist = add(env$selist, j  ,k, paste(s,e))
        }
        


        ## Recurse
        circBinSeg_inner(y=y,s=s,    e=s1, env=env)
        circBinSeg_inner(y=y,s=s1+1, e=s2, env=env)
        circBinSeg_inner(y=y,s=s2+1, e=e,  env=env)

    
    }
}


##' Computes the CBS CUSUM (hump) statistic.  Note, we calculate this as an
##' upward hump by default.
##' @param s starting index.
##' @param s1 first breakpoint index.
##' @param s2 second jbreakpoint index.
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
cusumm <- function(y,s=1,s1,s2,e=length(y), n=length(y), upward=TRUE, unsigned = FALSE, contrast.vec=FALSE){


    ## Basic checks
    if(s2>=e) stop("s2 must be strictly smaller than e!")
    if(s1>=s2) stop("s2 must be strictly larger than s1!")
    if(s>s1) stop("s1 must be larger than or equal to s!")
    if(unsigned & is.null(y)) stop("Can't produce an unsigned linear contrast or contrast vector without a y")
    if(is.null(y)){ if(is.null(n)) stop("Must provide one of y or n.")}
    if(is.null(n)) n = length(y)

    ## Calculate things
    ind.middle = (s1+1):s2
    ind.outer = c(s:s1, (s2+1):e)
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
    ssmat <- t(combn(s:(e-1),2))
    cusumms <- apply(ssmat, 1, function(ss){cusumm(y=y, s=s, s1=ss[1], s2=ss[2], e=e, unsigned = FALSE)})
    ind.max <- which.max(abs(cusumms))
    ss = ssmat[ind.max,]
    return(s1 = ss[1], s2 = ss[2], maxcusumm = cusumms[ind.max])
}


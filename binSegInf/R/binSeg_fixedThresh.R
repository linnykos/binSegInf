##' Main function for binary segmentation for fixed threshold. This is actually
##' a wrapper for binary segmentation with fixed threshold. It creates an
##' environment and creates the variables there, then runs
##' binseg.by.thresh.inner() all in this environment, and returns the relevant
##' guy
##' @param s Starting index, in the vector-valued data. Must be an integer
##'     larger than or equal to 1, and strictly smaller than \code{e}.
##' @param e Ending index, in the vector-valued data. Must be an integer smaller
##'     than or equal to \code{n}, and strictly larger than \code{s}.
##' @param thresh Threshold for \deqn{\tilde{X}}. This serves as a stopping rule
##'     for the recursion.
##' @param y The original data.
##' @param verbose Set to true if you'd like to see algorithm progression.
##' @param return.env Set to true if you'd like to get the environment
##'     containing the algorithm output, instead of a list.
##' @import pryr
binSeg_fixedThresh = function(y, thresh, s=1, e=length(y), verbose=FALSE, return.env=FALSE, numIntervals=NULL){

    n = length(y)
    
    ## Create new environment |env|
    env = new.env()
    env$slist = env$elist = env$selist = env$blist = env$Blist = env$zlist = env$Zlist =
        ## matrix(NA, nrow = n, ncol = 2 ^(n/2))
        ## matrix(NA,nrow=n,ncol=3)
        cplist(n)
    env$y = y

    ## Run binary segmentation on |env|
    binseg.by.thresh.inner(y, thresh, s, e, 0, 1, verbose, env=env)

    ## Gather output from |env| and return it.
    bs.output = list(slist = trim(env$slist),
                     elist = trim(env$elist),
                    selist = trim(env$selist),
                     blist = trim(env$blist),
                     Blist = trim(env$Blist),
                     ## zlist = trim(env$zlist),
                     Zlist = trim(env$Zlist),
                     cp = trim.vec(env$blist$mat[,3]),
                     cp.sign = trim.vec(env$zlist$mat[,3]),
                     y = y,
                     thresh = thresh)
    
    obj = structure(list(bs.output = bs.output,
                         y = y,
                         thresh = thresh,
                         cp = trim.vec(env$blist$mat[,3]),
                         cp.sign = trim.vec(env$zlist$mat[,3])) ,
                    class = "bsFt")
    
    if(return.env){ return(env) } else{ return(obj) }

}

##' Checks of object is of class "bsFt", produced by binSeg_fixedThresh().
##' @param obj bsFt object
is_valid.bsFt <- function(obj){
    if(!(all(names(obj) %in% c("bs.output",
                               "y",
                               "thresh",
                               "cp",
                               "cp.sign")))) stop("obj must contain certain elements!")
    TRUE
}



##' Inner function for binary segmentation with fixed threshold. The wrapper
##' \code{binseg.by.thresh()} is intended to be used by user. Note, when
##' \code{thresh} is set to zero, then this can be used to collect the
##' unbalanced haar wavelet basis.
##' @param s Starting index, in the vector-valued data. Must be an integer
##'     larger than or equal to 1, and strictly smaller than \code{e}.
##' @param e Ending index, in the vector-valued data. Must be an integer smaller
##'     than or equal to \code{n}, and strictly larger than \code{s}.
##' @param j The depth of the recursion on hand.
##' @param k The indexing of the node location, from left to right, in the
##'     \emph{complete} binary tree.
##' @param thresh Threshold for \deqn{\tilde{X}}. This serves as a stopping rule
##'     for the recursion.
##' @param y The original data.
##' @param n The length of the data \code{y}.

binseg.by.thresh.inner <- function(y, thresh, s=1, e=length(y), j=0, k=1, verbose=F, env=NULL){

    n = length(y)

    ## If segment is 2-lengthed, terminate
    if(e-s<1){
        env$slist = add(env$slist,j,k,s)
        env$elist = add(env$elist,j,k,e)
        env$selist = add(env$selist, j  ,k, paste(s,e))
        return() 
    ## Otherwise, calculate CUSUMs
    } else {
        all.bs = (s:(e-1))
        all.cusums = sapply(all.bs, function(b) cusum(s=s, b=b, e=e, y=y))
        names(all.cusums) = all.bs
        all.cusums = c(rep(NA,s-1), all.cusums, rep(NA,n-s))
        sn = sign(all.cusums) 
        
        ## Obtain breakpoint and its sign
        b = which.max(abs(all.cusums))
        z = sn[b]
        if(verbose) cat("the (potential) changepoint", b, "is being considered, between s and e:",s,e, fill=T)
        
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
        binseg.by.thresh.inner(y, thresh, s, b, j+1, 2*k-1, verbose, env=env)
        binseg.by.thresh.inner(y, thresh, b+1, e, j+1, 2*k, verbose, env=env)
    }
}



##' Print function for bsFt class
print.bsFt <- function(obj){
    ## if(obj$last.row==0){ print("Empty cplist object!")
    ## } else{ print(cplist$mat[1:cplist$last.row,])}
    
    cat("Changepoint set is ", obj$cp*obj$cp.sign,fill=TRUE)
}

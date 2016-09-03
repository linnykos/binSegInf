#' The main binary segmentation function. When \code{thresh} is set to zero,
#' then this can be used to collect the unbalanced haar wavelet basis.
#' @param s Starting index, in the vector-valued data. Must be an integer larger
#'     than or equal to 1, and strictly smaller than \code{e}.
#' @param e Ending index, in the vector-valued data. Must be an integer smaller
#'     than or equal to \code{n}, and strictly larger than \code{s}.
#' @param j The depth of the recursion on hand.
#' @param k The indexing of the node location, from left to right, in the
#'     \emph{complete} binary tree.
#' @param thresh Threshold for \deqn{\tilde{X}}. This serves as a stopping rule
#'     for the recursion.
#' @param y The original data.
#' @param n The length of the data \code{y}.
#' 
#' @examples
#' n = 60
#' sd = .5
#' mn = c(rep(-3, n/4), rep(2, n/4), rep(-1, n/4), rep(1, n/4))
#' y = mn + rnorm(n,0,sd)
#' blist = zlist = matrix(NA, nrow = n, ncol = 2^8)e
#' thresh = .5
#' binseg(s = 1, e = n, j = 0, k = 1, thresh = thresh, y = y, n = n)

binseg = function(s, e, j, k, thresh, y, n){
    cat(s,e,j,k,fill=T)
    if(e-s<1){
       ## cat(s,e,fill=T)
       return() 
    } else {
        leftindslist = lapply((s+1):e,
                              function(rightind){s:(rightind-1)}) 
        rightindslist = lapply((s+1):e,
                               function(leftind){leftind:e}) 

        ## Calculate all differences and their signs
        df = unlist(Map(function(leftind,rightind){ mean(y[rightind]) - mean(y[leftind]) },
                            leftindslist,
                            rightindslist))

        all.bs = (s:(e-1))
        df = sapply(all.bs, function(b) sqrt.mn.diff(s=s, b=b, e=e, n=n, y=y, contrast.vec=FALSE, right.to.left=TRUE))
        df = c(rep(NA,s-1),df, rep(NA,n-s))
        sn = sign(df) 

        ## Obtain breakpoint and its sign
        b = which.max(abs(df))
        z = sn[b]

        ## Check threshold exceedance, then store
        if(abs(df[b]) < thresh){
            return()
        } else {
            cat(df[b])
            blist[j+1,k] <<- b
            zlist[j+1,k] <<- z
        }
                    
        ## Recurse
        binseg(s, b, j+1, 2*k-1, thresh, y, n)
        binseg(b+1, e, j+1, 2*k, thresh, y, n)
    }
}



#' Helper function to get right-to-left mean difference, or the contrast vector
sqrt.mn.diff = function(s, b, e, n, y = NA, contrast.vec = FALSE,
                        right.to.left = TRUE){
    v = rep(0,n)
    v[s:b] = - sqrt(1/(b+1-s) - 1/(e-s+1))
    v[(b+1):e]  = + sqrt(1/(e-b) - 1/(e-s+1))
    if(!right.to.left) v = -v
    if(contrast.vec) return(v) else return(sum(v*y))
}




#' Calculates the unbalanced Haar basis for subvector of \code{y} or the vector
#' \code{a} in linear inequalities \( a y \ge 0 \).
#' @param type \code{basis} to return the unbalanced haar basis, and \code{ineq}
#'     to return linear inequality vector.
#' @examples
#' y = c(rnorm(10,0,1), rnorm(10,4,1))
#' 
#' # Calculate Haar basis, s:e vs (b+1):e
#' mybasis = haarbasis(s = 0, b = 10, e = 20, n = 20, y = y, type = "basis")
#' 
#' # Calculate linear inequality vectors
#' myineqs = haarbasis(s = 0, b = 10, e = 20, n = 20, y = y, type = "ineq")
#' 
haarbasis = function(s, b, e, n, y, type=c("basis", "ineq")){
    type = match.arg(type)
    
    if(!(s <= b & b <= e)){
        stop("s<=b<=e is not true")
    }

    ## Make contrast vector
    v = sqrt.mn.diff(s=s, b=b, e=e, n=n, contrast.vec=T, right.to.left=F)

    ## Make inequalities for halfspaces.
    V = matrix(0, nrow = n, ncol = n)
    other.bs = (s:(e-1))
    other.bs = other.bs[other.bs != b]
    ii = 0

    if(length(other.bs) == 0){
        V = V[c(),]
    } else {
        for(other.b in other.bs){
            cat(s,b,other.b,e,fill=T)
            ii = ii+1
            v.this = v.other = rep(NA,n)

            # Absolute mean difference at any other break.
            v.other = sqrt.mn.diff(s=s, b=other.b, e=e, n=n, contrast.vec=T,
                                   right.to.left=T)
            z.other = sign(sum(v.other * y))
            v.other = v.other * z.other

            # Absolute sqrt mean difference at maximizing break.
            v.this = sqrt.mn.diff(s=s, b=b, e=e, n=n, contrast.vec=T,
                                  right.to.left=T)
            z.this = sign(sum(v.this * y))
            v.this = v.this * z.this

            ## Quick sanity check of b
            stopifnot(sum(v.this  * y) > 0)
            stopifnot(sum(v.other * y) > 0)
            stopifnot(sum(v.this*y) > sum(v.other*y))

            # Subtract the two vectors, and store
            V[ii,] = v.other - v.this
        }
        V = V[1:ii,,drop=F]
    }

    return( (if(type=="basis") v else V))
}




#' Helper to get index of closest element of \code{val} out of vector \code{allval}
get.closest = function(val, allval){
    dists = (allval-val)
    pos.dists = neg.dists = dists
    pos.dists[dists<0] = Inf
    neg.dists[dists>0] = Inf 
    s = allval[which.min(abs(neg.dists))]
    e = allval[which.min(abs(pos.dists))]
    return(c(s+1,e))
}




#' Function to check if the columns in \code{basislist} are orthonormal.
#' @param tol Numerical allowance; letting inner products to be up to size
#'     \code{tol}.
check.orth.basis = function(basislist, tol = 1E-10){
    n = length(basislist)
    inner.products = matrix(NA,nrow=n,ncol=n) 
    for(i in 1:length(basislist)){
        for(j in 1:length(basislist)){
            cat(i,j,fill=T)
            if(i>j){
                inner.products[i,j] = (sum(basislist[[i]]*basislist[[j]]))
                print(inner.products[i,j])
            }
        }
    }
    print(inner.products)
    return(!any(as.numeric(inner.products)>tol, na.rm=T)) 
}

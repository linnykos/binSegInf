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
#' ## blist are the breaks, and zlist are the associated signs.
#' blist = zlist = matrix(NA, nrow = n, ncol = 2^8)
#' thresh = .5
#' binseg(s = 1, e = n, j = 0, k = 1, thresh = thresh, y = y, n = n)
binseg = function(s, e, j, k, thresh, y, n){
    cat("s,e,j,k",fill=T)
    cat(s,e,j,k,fill=T)
    if(e-s<1){
        cat("terminated because e-s<1",fill=T)
       slist[j,k] <<- s
       elist[j,k] <<- e
       return() 
    } else {
        all.bs = (s:(e-1))
        df = sapply(all.bs, function(b) sqrt.mn.diff(s=s, b=b, e=e, n=n, y=y, contrast.vec=FALSE, right.to.left=TRUE))
        df = c(rep(NA,s-1),df, rep(NA,n-s))
        sn = sign(df) 

        ## Obtain breakpoint and its sign
        b = which.max(abs(df))
        z = sn[b]
        cat("b is", b, fill=T)

        ## Check threshold exceedance, then store
        if(abs(df[b]) < thresh){
            Blist[j+1,k] <<- b
            Zlist[j+1,k] <<- z
            slist[j,k] <<- s
            elist[j,k] <<- e
            cat("terminated because biggest gap was",abs(df[b]),fill=T)
            return()
        } else { 
            cat(df[b], fill=T)
            blist[j+1,k] <<- b
            Blist[j+1,k] <<- b
            zlist[j+1,k] <<- z
            Zlist[j+1,k] <<- z
            slist[j,k] <<- s
            elist[j,k] <<- e
        }
                    
        ## Recurse
        binseg(s, b, j+1, 2*k-1, thresh, y, n)
        binseg(b+1, e, j+1, 2*k, thresh, y, n)
    }
}



#' Calculates the halfspace vector for the maximizing breakpoint and all the signs.
#' @examples
#' y = c(rnorm(10,0,1), rnorm(10,4,1))
#' # Calculate linear inequality vectors
#' myineqs = basis.hs(s = 0, b = 10, e = 20, n = 20, y = y, type = "ineq")
halfspaces = function(s, b, e, thresh, n, y, terminal=F){
    
    if(!(s <= b & b <= e)){
        stop("s<=b<=e is not true")
    }

    V = matrix(NA, nrow = n^2, ncol = n)
    u = rep(NA,n)
    other.bs = (s:(e-1))
    other.bs = other.bs[other.bs != b]
    ii = 0

    ## Sqr mean difference of (s,b,e)
    v.this = sqrt.mn.diff(s=s, b=b, e=e, n=n, contrast.vec=T,
                          right.to.left=T)
    z.this = sign(sum(v.this * y))
    vz.this = v.this * z.this

    ## Characterizing this break's sign.
    ii = ii+1
    V[ii,] = vz.this
    u[ii] = 0

    ## Characterizing Gap size exceedance
    if(!terminal){
        ii = ii+1
        V[ii,] = vz.this## (if(!terminal) vz.this else -vz.this)
        u[ii] = thresh ##(if(!terminal) thresh else -thresh)
    }

    if(length(other.bs) == 0){
        V = V[c(),]
    } else {
        for(other.b in other.bs){
           ## Sqrt mean difference of (s,other.b,e)
            v.other = sqrt.mn.diff(s=s, b=other.b, e=e, n=n, contrast.vec=T,
                                   right.to.left=T)
            z.other = sign(sum(v.other * y))
            vz.other = v.other * z.other

            ## Characterizing other breaks' signs
            ii = ii+1
            V[ii,] = vz.other
            u[ii] = 0

            ## Characterizing maximizer of |sqrt mean difference|
            if(!terminal){
                ii = ii+1
                V[ii,] = -vz.other + vz.this
                u[ii] = 0
            }
            
            ## Characterizing failure to exceed gap size (for terminal nodes)
            if(terminal){
                ii = ii + 1
                V[ii,] = -vz.other
                u[ii] = -thresh
            }
                
            ## Quick sanity check of b
            stopifnot(sum(vz.this  * y) > 0)
            stopifnot(sum(vz.other * y) > 0)
            stopifnot(sum(vz.this*y) > sum(vz.other*y))
            if(!terminal){ stopifnot(sum(vz.this*y) > thresh)}
            if(terminal){  stopifnot(sum(-vz.other*y) > -thresh)}
        }
        V = V[1:ii,,drop=F]
        u = u[1:ii]
    }
    return(list(V=V,u=u))
}


#' Helper to see if V is big enough to store ii'th row. 
#' @examples
#' if(!big.enough(V,ii)){V = rbind(V, matrix(NA,nrow=nrow(V),ncol=n)}
big.enough = function(V,ii){
    nrow(V)
    ## Not written yet.
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
haarbasis = function(s, b, e, n, y, type=c("basis", "ineq"),){
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



#' Helper to collapse matrix (with NAs) to a vector of unique elements.
    collapse.prev = function(mat){
        collapsed = as.numeric(mat)
        collapsed = collapsed[!is.na(collapsed)]
        collapsed = unique(collapsed)
        return(collapsed)
    }


#' Function to trim a matrix from the right and bottom, ridding of all-NA rows/columns.
trim = function(mat, type = c("rowcol","row")){
    type = match.arg(type)
    if(is.null(dim(mat))){ mat = mat[1:max(which(!is.na(mat)))]; return(mat)}
    last.j = max(which(!(apply(mat,1,function(myrow) return(all(is.na(myrow)))))))
    mat = mat[1:last.j,,drop=F]
    if(type=="rowcol"){
        last.j = max(which(!(apply(mat,2,function(mycol) return(all(is.na(mycol)))))))
        mat = mat[,1:last.j,drop=F]
    }
    return(mat)
}

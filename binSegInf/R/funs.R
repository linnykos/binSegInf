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
#' ##binseg(s = 1, e = n, j = 0, k = 1, thresh = thresh, y = y, n = n) ## check() doesn't like this; not sure why.
binseg = function(s, e, j, k, thresh, y, n, verbose=F){
    if(verbose){
        cat("s,e,j,k",fill=T)
        cat(s,e,j,k,fill=T)
    }
    if(e-s<1){
       if(verbose) cat("terminated because e-s<1", fill=T)
       slist[j,k] <<- s
       elist[j,k] <<- e
       return() 
    } else {
        all.bs = (s:(e-1))
        df = sapply(all.bs, function(b) cusum(s=s, b=b, e=e, y=y))
        df = c(rep(NA,s-1),df, rep(NA,n-s))
        sn = sign(df) 

        ## Obtain breakpoint and its sign
        b = which.max(abs(df))
        z = sn[b]
        if(verbose) cat("b is", b, fill=T)

        ## Check threshold exceedance, then store
        if(abs(df[b]) < thresh){
            Blist[j+1,k] <<- b
            Zlist[j+1,k] <<- z
            slist[j,k] <<- s
            elist[j,k] <<- e
            if(verbose) cat("terminated because biggest gap was",abs(df[b]),fill=T)
            return()
        } else { 
            if(verbose) cat(df[b], fill=T)
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




#' Calculates the halfspace vectors for the maximizing breakpoint and all the
#' signs.
#' @param is.terminal.node T/F for whether the node is one where the threshhold
#'     is not breached.
#' @examples
#' y = c(rnorm(10,0,1), rnorm(10,4,1))
#' # Calculate linear inequality vectors
#' ##myineqs = halfspaces(s = 0, b = 10, e = 20, n = 20, y = y, type = "ineq")
halfspaces = function(s, b, e, thresh, n, y, is.terminal.node=F , verbose=F){
    
    if(verbose) cat("s,b,e are", s,b,e, fill=T)
    if(!(s <= b & b <= e)){
        stop("s<=b<=e is not true")
    }

    V = matrix(NA, nrow = n^2, ncol = n)
    u = rep(NA,n)
    other.bs = (s:(e-1))
    other.bs = other.bs[other.bs != b]
    ii = 0

    ## CUSUM comparison (s,b,e)
    v.this = cusum(s=s, b=b, e=e, n=n, contrast.vec=TRUE, right.to.left=TRUE)
    z.this = sign(sum(v.this * y))
    vz.this = v.this * z.this

    ## Characterizing this break's sign.
    ii = ii+1
    V[ii,] = vz.this
    u[ii] = 0

    ## Characterizing Gap size exceedance
    if(!is.terminal.node){
        ii = ii+1
        V[ii,] = vz.this
        u[ii] = thresh
    }

    if(length(other.bs) == 0){
        V = V[c(),]
    } else {
        for(other.b in other.bs){
           ## Sqrt mean difference of (s,other.b,e)
            v.other = cusum(s=s, b=other.b, e=e, n=n, contrast.vec=T,
                            right.to.left=T)
            z.other = sign(sum(v.other * y))
            vz.other = v.other * z.other

            ## Characterizing other breaks' signs
            stopifnot(sum(vz.this  * y) > 0)
            stopifnot(sum(vz.other * y) > 0)
            ii = ii+1
            V[ii,] = vz.other
            u[ii] = 0

            ## Characterizing maximizer of |sqrt mean difference|
            if(!is.terminal.node){
                stopifnot(sum(vz.this*y) > sum(vz.other*y))
                ii = ii+1
                V[ii,] = -vz.other + vz.this
                u[ii] = 0
            }
            
            ## Characterizing failure to exceed gap size (for terminal nodes)
            if(is.terminal.node){
                stopifnot(sum(-vz.this*y) > -thresh)
                ii = ii + 1
                V[ii,] = -vz.other
                u[ii] = -thresh
            }
        }

        V = V[1:ii,,drop=F]
        u = u[1:ii]
    }
    return(list(V=V,u=u))
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

#' Function to collect polyhedrons given some output from the binary
#' segmentation function. Takes as input list contatining *list objects. 
get.polyhedron = function(binseg.results, thresh, verbose = F) {

    ## Extracting things
    Blist = binseg.results$Blist
    blist = binseg.results$blist
    slist = binseg.results$slist
    elist = binseg.results$elist
    zlist = binseg.results$zlist
    Zlist = binseg.results$Zlist
    y = binseg.results$y
    thresh = binseg.results$thresh

    ## Initialize G and u
    ii = 1
    G = matrix(NA,nrow = nrow(trim(Blist))*n, ncol = n)
    u = rep(NA, nrow(trim(Blist))*n)
    nrow.G = 0

    ## For each node in Blist, collect halfspaces.
    for(my.j in 1:nrow(trim(Blist))){
        if(verbose) cat('\r', 'step', my.j)
        bsublist = Blist[my.j,]
        for(b in bsublist[!is.na(bsublist)]){
            if(my.j==1){k=1
                se = c(1,n)
            } else {
                k = which(b==bsublist)
                se = c(slist[my.j-1,
                             k], elist[my.j-1,k]) 
            }
            
            ## Collect halfspaces
            my.halfspaces = halfspaces(s = se[1],
                                       e = se[2],
                                       b = b,
                                       thresh = thresh,
                                       n = n,
                                       y = y,
                                       is.terminal.node = is.na(blist[my.j,k]))
            newrows = my.halfspaces[["V"]]
            newconst = my.halfspaces[["u"]]
            
            ## Move on if no comparisons were made i.e. lengths to left=2, right=1
            if(dim(newrows)[1]==0) next
            
            ## Add to Gamma matrix
            newrowinds = nrow.G + c(1:nrow(newrows))
            if(any(newrowinds > nrow(G) )){
                G = rbind(G, matrix(NA,nrow=nrow(G),ncol=n))
                u = c(u, rep(NA,length(u)))
            }
            G[newrowinds,] = newrows
            u[newrowinds] = newconst 
            
            ## Updates for loop
            nrow.G = nrow.G + nrow(newrows)
            ii = ii + 1
        }
    }
    if(verbose) cat(fill=T)

    ## Trim and return
    G = trim(G,"row")
    u = trim(u)
    return(list(G=G, u=u))
}

##' Get all cusums, given S and E
getcusums = function(s,e,y,returnmax = TRUE){
   cusums =  sapply(s:(e-1), function(b){cusum(s,b,e,y)})
   if(returnmax){return(list(bmax = which.max(cusums)+s-1, cusum = max(cusums)))}
   return(cusums)
}


##' Function to carry out fixed-number-of-steps binary segmentation.
##' @param y numeric vector, data
##' @param k desired number of changepoints
##' @import Matrix

binseg.by.thresh = function(y,k){
    ## Example
    set.seed(0) 
    y = c(rnorm(10,0,.5),rnorm(10,3,.5), rnorm(10,5,.5))
    numsteps = 10
    n=30
    
    ## initialize thin
    ## A=T=S=E=list()
    B = rep(NULL,n) 

    Bcurr=Scurr=Ecurr=Matrix(0,ncol=2^numsteps, nrow = numsteps,sparse=TRUE)
    A=T=S=E=Tcurr=Acurr=lapply(1:n,function(i) c())
    jk  = list()

    zetas = c()
    Scurr[1,1] = 1
    Ecurr[1,1] = length(y)
    Tcurr[[1]] = c(1,1)
    Acurr[[1]] = c()
    
    for(mystep in 1:numsteps){
        print(mystep)

        ## Get all candidate changepoints
        curr.max=-Inf
        for(ii in 1:sum(!sapply(Tcurr, is.null))){ 
            j = Tcurr[[ii]][1]
            k = Tcurr[[ii]][2]
            cusums = getcusums(s = Scurr[j,k],
                               e = Ecurr[j,k],
                               y = y,
                               returnmax = TRUE)
            Bcurr[j, k] = cusums$bmax
            breaking.cusum = cusums$cusum
            
            ## Keep running maximum
            if(curr.max <= breaking.cusum){
                curr.which.max = c(j,k) 
                curr.max = breaking.cusum
            }
        }
        
        ## Record knot as CUSUM maximizer
        zetas[mystep] = curr.max
        jmax = curr.which.max[1]
        kmax = curr.which.max[2]
        
        ## Update terminal and active node set
        which.duplicate = which(sapply(Tcurr, function(myjk){all.equal(myjk, c(jmax, kmax))==TRUE}))
        Tcurr[[which.duplicate]] <- c(jmax + 1, 2*kmax - 1)
        Tcurr[[mystep+1]] <- c(jmax + 1, 2*kmax)
        Acurr[[mystep+1]] <- c(jmax,kmax)
        
        ## Update Scurr and Ecurr for the /new/ nodes
        Scurr[jmax+1,2*kmax-1] = Scurr[jmax,kmax]
        Ecurr[jmax+1,2*kmax-1] = Bcurr[jmax,kmax]
        
        Scurr[jmax+1,2*kmax] = Bcurr[jmax,kmax]+1
        Ecurr[jmax+1,2*kmax] = Ecurr[jmax,kmax]

        ## Take snapshot
        S[[mystep]] = trim(Scurr)
        E[[mystep]] = trim(Ecurr)
        A[[mystep]] = trimlist(Acurr)
        T[[mystep]] = trimlist(Tcurr)
        B[mystep] = Bcurr[jmax,kmax]
        jk[[mystep]] = c(jmax,kmax) 
    }

    ## trim and return things
    return(list(S = trimlist(S),
                E = trimlist(E),
                A = trimlist(A),
                T = trimlist(T),
                B = trimlist(B)))
}    


##' Computes the CUSUM (cumulative sum) statistic. Note, we calculate this as
##' the right-to-left difference, by default.
##' @param s starting index.
##' @param b breakpoint index.
##' @param e end index.
##' @param y data.

cusum = function(s,b,e,y, right.to.left = TRUE){

    ## Form temporary quantities
    n = e-b+1
    n1 = b+1-s
    n2 = e-b
    if(n==1) stop("n cannot be 1!")
    if(b>=e) stop("b must be strictly smaller than e!")
    if(s>b) stop("b must be larger than or equal to!")

    ## Form contrast
    v[s:b] = -1/n1 
    v[(b+1):e]  = 1/n2
    v = v * sqrt(1/((1/n1)+(1/n2)))
    if(!right.to.left) v = -v

    ## Return the right thing
    if(contrast.vec) return(v)
    return(sum(v*y))
}



#' Function to obtain contrast
make.v = function(test.b, bs.output){

    ## Extract values
    G = bs.output$G
    u = bs.output$u
    y = bs.output$y
    blist = bs.output$blist
    zlist = bs.output$zlist

    ## Basic checks
    stopifnot(test.b %in% sort(collapse(trim(bs.output$blist))))
    stopifnot(all(!is.na(collapse(G))))
    stopifnot(all(!is.na(collapse(u))))
    
    ## Make contrast vector
    n = length(bs.output$y)
    v = rep(0,n)
    z = zlist[which(test.b == blist, arr.ind = T)]
    ends = c(0,sort(blist),n)
    ind = which(test.b == ends)
    my.se = ends[c(ind-1, ind+1)]
    left.b = (my.se[1]+1):(test.b)
    right.b = (test.b+1):(my.se[2])
    v[left.b] = -1/length(left.b)
    v[right.b] = 1/length(right.b)
    v = v*z

    ## Only for one-sided tests
    stopifnot(v%*%y>0)
    
    return(v)
}



#' Calculates the unbalanced Haar basis for subvector of \code{y} or the vector
#' \code{a} in linear inequalities a*y>=0
#' @param type \code{basis} to return the unbalanced haar basis, and \code{ineq}
#'     to return linear inequality vector.
#' @examples
#' y = c(rnorm(10,0,1), rnorm(10,4,1))
#' 
#' # Calculate Haar basis, s:e vs (b+1):e
#' ##mybasis = haarbasis(s = 0, b = 10, e = 20, n = 20, y = y, type = "basis")
#' 
#' # Calculate linear inequality vectors
#' ##myineqs = haarbasis(s = 0, b = 10, e = 20, n = 20, y = y, type = "ineq")
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
collapse = function(mat){
    collapsed = as.numeric(mat)
    collapsed = collapsed[!is.na(collapsed)]
    collapsed = unique(collapsed)
    return(collapsed)
}


##' Function to trim a matrix from the right and bottom, ridding of all-NA rows/columns.
##' Returns NULL if mat is all NA's.
trim = function(mat, type = c("rowcol","row")){
    type = match.arg(type)

    ## If all NA matrix, return NULL.
    if(all(is.na(as.numeric(mat)))) return(NULL)
    
    if(is.null(dim(mat))){ mat = mat[1:max(which(!is.na(mat)))]; return(mat)}
    last.j = max(which(!(apply(mat,1,function(myrow) return(all(is.na(myrow)))))))
    mat = mat[1:last.j,,drop=F]
    if(type=="rowcol"){
        last.j = max(which(!(apply(mat,2,function(mycol) return(all(is.na(mycol)))))))
        mat = mat[,1:last.j,drop=F]
    }
    return(mat)
}

trimlist = function(mylist){
    return(mylist[1:(max(which(!sapply(mylist, is.null))))])
}


#' Calculates a t-statistic for E(v2)-E(v1) given vectors v1 and v2.
t.statistic = function(v1, v2){
    numer = (mean(v2)-mean(v1))
    denom = sqrt(var(v1)/length(v1) + var(v2)/length(v2))
    return(numer/denom) 
}

#' Conducts a permutation t-test, given two subvectors
#' Example doesn't work now, but here it is: EXAMPLES/perm.t.test.example.R
perm.t.test = function(vec1, vec2, nsim=1000){
    vec = c(vec1,vec2)
    n1 = length(vec1)
    n2 = length(vec2)
    n = n1 + n2
    null.tstats = replicate(nsim,{
        newvec = vec[sample(1:n,n,replace=FALSE)]
        newvec1 = newvec[1:n1]
        newvec2 = newvec[(n1+1):(n2)]
        t.statistic(newvec1,newvec2)})
    pval = 1-ecdf(null.tstats)(t.statistic(vec1,vec2))
    cat('Ran', nsim, 'permuatations for the permutation t-test.', fill=TRUE)
    return(pval)
}


##' Gets underlying means for changepoints in y.  Wrote this because plot.sbs
##' function doesn't work properly.
##' Example doesn't work now, but here it is: EXAMPLES/get.mean.example.R
get.means = function(y, changepoints, ...){
    cps = c(0, sort(changepoints), n)
    mns = rep(NA, length(y))
    for(ii in 1:(length(cps)-1) ){
        inds = (cps[ii]+1):(cps[ii+1])
        mns[inds] = mean(y[inds])
    }
    if(any(is.na(mns))) stop("NAs in mns!")
    return(mns)
}



# Makeshift replacement of all usages of pval.fl1d to poly.pval; 
pval.fl1d <- function(y, G, dik, sigma, approx=T, threshold=T, approxtype = c("gsell","rob"), u = rep(0,nrow(G))){
  return(poly.pval(y, G, u, dik, sigma, bits=NULL)$pv)
}

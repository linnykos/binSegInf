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
getcusums = function(s,e,y){

    if(s<=0 | e<= 0) stop("must enter valid e,s >=1 ")

    ## Get all cusum
    cusums =  sapply(s:(e-1), function(b){cusum(s,b,e,y)})
    names(cusums) = paste("b=",s:(e-1))
    contrasts =  t(sapply(s:(e-1), function(b){cusum(s,b,e,y, contrast.vec=TRUE)}))

    ## Get signs
    signs = sign(cusums)
    ## abs.cusums = abs(cusums)
    abs.cusums = signs*cusums

    return(list(bmax = which.max(abs.cusums)+s-1, bmax.cusums = which.max(abs.cusums),
                cusum = max(abs.cusums), allcusums = cusums, contrasts = contrasts,
                signs=signs))
}


##' Computes the CUSUM (cumulative sum) statistic. Note, we calculate this as
##' the right-to-left difference, by default.
##' @param s starting index.
##' @param b breakpoint index.
##' @param e end index.
##' @param y data.
##' @param right.to.left Whether you want right-to-left difference in the cusum calculation. Defaults to TRUE.
##' @param contrast.vec If TRUE, then the contrast vector v for cusum=v'y is returned.

cusum = function(s,b,e,y, right.to.left = TRUE, contrast.vec = FALSE){

    ## Form temporary quantities
    n = e-b+1
    n1 = b+1-s
    n2 = e-b
    if(n==1) stop("n cannot be 1!")
    if(b>=e) stop("b must be strictly smaller than e!")
    if(s>b) stop("b must be larger than or equal to!")

    ## Form contrast
    v = rep(0,length(y))
    v[s:b] = -1/n1 
    v[(b+1):e]  = 1/n2
    v = v * sqrt(1/((1/n1)+(1/n2)))
    if(!right.to.left) v = -v

    ## Return the right thing
    if(contrast.vec) return(v)
    return(sum(v*y))
}




##' Function to obtain contrast, given the output of binseg(), the SBS algorithm
##' with fixed threshold.
##' @param test.b is the location that we want to test.
##' @param bs.output list that contains G,u,y,blist,zlist. Manually bundled by
##'     user.
make.v.fixed.thresh = function(test.b, bs.output){

    ## Extract values
    y = bs.output$y  ## 
    blist = bs.output$blist
    zlist = bs.output$zlist

    ## Basic checks
    stopifnot(test.b %in% sort(collapse(trim(blist))))
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


##' Function to obtain contrast, given the output of binseg(), the SBS algorithm
##' with fixed threshold.
##' @param test.b is the location that we want to test.
##' @param B vector of breakpoints to be considered
##' @param Z signs of B, as 
##' @param n length of contrast 
make.v = function(test.b,B,Z,n){
        v = rep(0,n)
        z = Z[which(test.b == B, arr.ind = T)]
        ends = c(0,sort(B),n)
        ind = which(test.b == ends)
        my.se = ends[c(ind-1, ind+1)]
        left.b = (my.se[1]+1):(test.b)
        right.b = (test.b+1):(my.se[2])
        v[left.b] = -1/length(left.b)
        v[right.b] = 1/length(right.b)
        v = v*z
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
trim.mat = function(mat, type = c("rowcol","row")){
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

##' Trims a list by deleting the last consecutive elements that are NULL.
##' @param mylist Some list.
trim.list = function(mylist){
    return(mylist[1:(max(which(!sapply(mylist, is.null))))])
}

##' Trims a list by deleting the last consecutive elements that are NULL.
##' @param mylist Some list.
trim.vec = function(myvec){
    return(myvec[1:(max(which(!sapply(myvec, is.na))))])
}


##' Function to trim matrices, lists or vectors.
trim = function(mything,...){
    class.of.my.thing = class(mything)
    if(class.of.my.thing == "list"){
        return(trim.list(mything,...))
    } else if (class.of.my.thing %in% c("matrix","dgCMatrix")){
        return(trim.mat(mything,...))
    } else if (class.of.my.thing %in% c("integer", "numeric")){
        return(trim.vec(mything,...))
    } else {
        stop(paste("trim() doesn't know how to trim things of class:", class.of.my.thing))
    }
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



##' Takes a contrast vector v and optionally B & Z, and plots it
##' @param v Numeric vector containing contrast vector.
plot.v = function(v, B=NULL, Z=NULL){
    ## Basic checks
    if(any(is.na(B))|any(is.na(Z))|length(Z)!=length(B)) stop("Check yo Z and B!")
    stopifnot(all(abs(v)<=1))
    stopifnot(all(Z%in%c(-1,1)))

    ## Make plot
    maxabs = max(abs(v))
    plot(NA,ylim=c(-1,+1)*maxabs*1.3, xlim = c(0,length(v)))
    lines(v, type='o', col = 'blue', pch = 16)
    if(!is.null(B))  abline(v=B, col='grey70', lty=2)
    if(!is.null(Z)) points(x=B, y=rep(maxabs,length(B)),
                           pch = sapply(Z,function(myz){if(myz==+1)"+"else"-"}))
}

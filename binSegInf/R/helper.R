##' Computes the CUSUM (cumulative sum) statistic.
##'
##' Note, we calculate this as the right-to-left difference, by default.
##' @param s starting index.
##' @param b breakpoint index.
##' @param e end index.
##' @param y data.
##' @param n length of data. Defaults to \code{length(y)}.
##' @param right.to.left Whether you want right-to-left difference in the cusum
##'     calculation. Defaults to TRUE.
##' @param contrast.vec If TRUE, then the contrast vector v for cusum=v'y is
##'     returned.
##' @param unsigned if TRUE, then returns the contrast vector that makes the
##'     contrast with |y| positive, or the absolute value of the contrast.
##' @export

cusum <- function(s,b,e,n=NULL, y=NULL, right.to.left = TRUE, contrast.vec = FALSE, unsigned = FALSE){

    ## Form temporary quantities
    nn = e-b+1
    n1 = b+1-s
    n2 = e-b
    if(nn==1) stop("n cannot be 1!")
    if(b>=e) stop("b must be strictly smaller than e!")
    if(s>b) stop("b must be larger than or equal to!")
    if(unsigned & is.null(y)) stop("Can't produce an unsigned linear contrast or contrast vector without a y")
    ## if(!is.null(y)) stopifnot(length(y)==n)
    if(is.null(y)){ if(is.null(n)) stop("Must provide one of y or n.")}
    if(is.null(n)) n = length(y)

    ## Form contrast
    v = rep(0,n)
    v[s:b] = -1/n1 
    v[(b+1):e]  = 1/n2
    v = v * sqrt(1/((1/n1)+(1/n2)))
    if(!right.to.left) v = -v
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


#' Calculates the halfspace vectors for the maximizing breakpoint and all the
#' signs, for fixed-threshold SBS.
#' 
#' @param is.terminal.node T/F for whether the node is one where the threshold
#'     is not breached.
#' 
#' @return A list of two objects \code{V} and \code{u}, for the halfspaces in
#'     represented as V'y>u
#' 
#' @examples 
#' y = c(rnorm(10,0,1), rnorm(10,4,1))
#' myineqs = halfspaces(s = 0, b = 10, e = 20, n = 20, y = y)
halfspaces = function(s, b, e, z, thresh, n, y, is.terminal.node=F , verbose=F){
    if(verbose) cat("s,b,e are", s,b,e, fill=T)
    if(!(s <= b & b <= e)){
        stop("s<=b<=e is not true")
    }
    
    ## Make empty things, initialize
    V = matrix(NA, nrow = n^2, ncol = n)
    u = rep(NA,n)
    other.bs = (s:(e-1))
    other.bs = other.bs[other.bs != b]
    ii = 0
    
    ## CUSUM comparison (s,b,e)
    v.this = cusum(s=s, b=b, e=e, y=y, contrast.vec=TRUE, right.to.left=TRUE)
    ## z.this = sign(sum(v.this * y))
    z.this = z
    vz.this = v.this * z.this
    
    ## 1. Characterizing this break's sign.
    ii = ii+1
    V[ii,] = vz.this
    u[ii] = 0
    
    ## 2. Characterizing threshold exceedance/nonexceedance
    if(!is.terminal.node){
        ii = ii+1
        V[ii,] = vz.this
        u[ii] = thresh
    } else {
        stopifnot(sum(-vz.this*y) > -thresh)
        ii = ii + 1
        V[ii,] = -vz.this
        u[ii] = -thresh
    }
    
    if(length(other.bs) == 0){
        V = V[c(),]
    } else {
        for(other.b in other.bs){
            
            ## Sqrt mean difference of (s,other.b,e)
            v.other = cusum(s=s, b=other.b, e=e, y=y, contrast.vec=T,
                            right.to.left=T)
            
            ## 3. Characterizing maximizer of |sqrt mean difference|
            if(!is.terminal.node){
                ## other cusum is larger than -|this cusum|
                ii = ii+1
                stopifnot(-sum(vz.this*y) < sum(v.other*y))
                V[ii,] = vz.this + v.other
                u[ii] = 0
                ## other cusum is smaller than +|this cusum|
                ii = ii+1
                stopifnot(sum(v.other*y) < sum(vz.this*y) )
                V[ii,] = vz.this - v.other
                u[ii] = 0
            }
        }
        
        V = V[1:ii,,drop=F]
        u = u[1:ii]
    }
    return(list(V=V,u=u))
}



##' From polyhedron and data vector and contrast, gets the probability that vtY
##' is in the polyhedron, conditional on PvperpY. i.e. the probability that Y
##' stays in the polyhedron, fixing n-1 dimensions of it.
##' @param y data
##' @param polyhedra polyhedra object produced form \code{polyhedra(wbs_object)}
##' @param sigma data noise (standard deviation)
##' @param nullcontrast the null value of \eqn{v^T\mu}, for \eqn{\mu = E(y)}.
##' @param v contrast vector
##'
##' @return list of two vectors: denominators and numerators, each named
##'     \code{denom} and \code{numer}.
partition_TG <- function(y, polyhedra, v, sigma, nullcontrast=0, bits=50){
    
    ## Basic checks
    stopifnot(length(v)==length(y))
    stopifnot(is_valid.polyhedra(polyhedra))
    
    ## From selectiveinference package (todo: make it succinct)
    G = polyhedra$gamma
    u = polyhedra$u
    mean=0
    
    z = sum(v*y)
    vv = sum(v^2)
    sd = sigma*sqrt(vv)
  
    rho = G %*% v / vv
    vec = (u - G %*% y + rho*z) / rho
    a = vlo = suppressWarnings(max(vec[rho>0]))
    b = vup = suppressWarnings(min(vec[rho<0]))
    

    ## Calculate a,b,z for TG = (F(b)-F(z))/(F(b)-F(a))
    z = Rmpfr::mpfr((z-mean)/sd, precBits=bits)
    a = Rmpfr::mpfr((a-mean)/sd, precBits=bits)
    b = Rmpfr::mpfr((b-mean)/sd, precBits=bits)

    ## Separately store and return num&denom of TG
    numer = as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z)))
    denom = as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(a)))
    return(list(denom=denom, numer=numer))
}



##' Function to plot qqlot of p-values
##' @param pp numeric vector of p-values.
##' @param main 
qqunif <- function(pp, main=NULL){
    qqplot(x=pp, y=seq(from=0,to=1,length=length(pp)))
    abline(0,1)
    if(!is.null(main)) title(main=main)
}



##' Get all cusums, given start point \code{s} and end point \code{e}
##' 
##' @param s Starting index, between \code{1} and \code{n}
##' @param e Ending index, between \code{1} and \code{n}
##' @param y \code{n}-lengthed data vector.
##' @param unsigned \code{TRUE} to return alsolute value.
##'
##' @return list of information about the cusum calculations and maximizers in
##'     this interval.
##' @example examples/getcusums-example.R
##' @export
getcusums <- function(s,e,y, unsigned=FALSE){

    if(s<=0 | e<= 0) stop("must enter valid e,s >=1 ")

    ## Get all cusum
    cusums =  sapply(s:(e-1), function(b){ cusum(s=s,b=b,e=e,y=y, unsigned=unsigned) })
    names(cusums) = paste("b=",s:(e-1))
    contrasts =  t(sapply(s:(e-1), function(b){cusum(s=s,b=b,e=e,y=y, contrast.vec=TRUE, unsigned = unsigned)}))

    ## Get signs
    signs = sign(cusums)
    abs.cusums = signs*cusums

    return(list(bmax = which.max(abs.cusums)+s-1,
                bmax.cusums = which.max(abs.cusums),
                inds = s:(e-1),
                cusum = max(abs.cusums),
                allcusums = cusums,
                contrasts = contrasts,
                signs=signs))
}


##' Either return maximizing breakpoint, or the maximum cusum
.get_max_b <- function(m,intervals,y,type=c("cusums","max.b", "max.z")){
    type = match.arg(type)
    s <- intervals$starts[[m]]
    e <- intervals$ends[[m]]
    max.b.before.adding.s <- which.max(getcusums(s=s,e=e,y=y,unsigned=TRUE)$allcusums)
    if(type=="cusums"){
        return(max(getcusums(s=s,e=e,y=y,unsigned=TRUE)$allcusums))
    } else if (type == "max.b"){
        ## max.b.before.adding.s <- which.max(getcusums(s=s,e=e,y=y,unsigned=TRUE)$allcusums)
        max.b <- max.b.before.adding.s + (s - 1)
        return(max.b)
    } else if (type == "max.z"){
        allcusums = getcusums(s=s,e=e,y=y,unsigned=FALSE)$allcusums
        max.zs <- sign(allcusums)
        max.z <- max.zs[max.b.before.adding.s]
        return(max.z)
    } else {
        stop("type is not defined")
    }
}


##' Function to crete 1d fused lasso regression matrix.
##' @param m data length
##' @return 1d Fused lasso regression matrix
dual1d_Dmat = function(m){
  D = matrix(0, nrow = m-1, ncol = m)
  for(ii in 1:(m-1)){
    D[ii,ii] = -1
    D[ii,ii+1] = 1
  }
  return(D)
}

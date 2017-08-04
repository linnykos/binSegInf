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
    if(unsigned) v = v * sign(sum(v*y)) ## This needs to change! Have the
                                        ## function require the sign of v*y
                                        ## beforehand, or check only if it isn't
                                        ## provided, or something!

    ## Return the right thing
    if(contrast.vec){
        return(v)
    } else if(!is.null(y)) {
        return(sum(v*y))
    } else {
        stop("Check your options for cusum() again!")
    }
}


#' Calculates the halfspace vectors for the maximiiming breakpoint and all the
#' signs, for fixed-threshold SBS.
#'
#' @param is.terminal.node T/F for whether the node is one where the threshold
#'     is not breached.
#'
#' @return A list of two objects \code{V} and \code{u}, for the halfspaces in
#'     represented as V'y>u

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
##' @param poly polyhedra object produced form \code{polyhedra(wbs_object)}
##' @param sigma data noise (standard deviation)
##' @param nullcontrast the null value of \eqn{v^T\mu}, for \eqn{\mu = E(y)}.
##' @param v contrast vector
##'
##' @return list of two vectors: denominators and numerators, each named
##'     \code{denom} and \code{numer}.
partition_TG <- function(y, poly, v, sigma, nullcontrast=0, bits=50, reduce){

    ## Basic checks
    stopifnot(length(v)==length(y))
    ## stopifnot(is_valid.polyhedra(poly))

    vy = sum(v*y)
    vv = sum(v^2)
    sd = sigma*sqrt(vv)

    ## Just in case |poly| doesn't contain |vup| and |vlo|, we manually form it.
    ## This is because in order to partition the TG statistic, we need to form
    ## these anyway.
    pvobj <- poly.pval2(y, poly, v, sigma)
    vup = pvobj$vup
    vlo = pvobj$vlo
    vy = max(min(vy, vup),vlo)

    ## Calculate a,b,z for TG = (F(b)-F(z))/(F(b)-F(a))
    z = Rmpfr::mpfr(vy/sd, precBits=bits)
    a = Rmpfr::mpfr(vlo/sd, precBits=bits)
    b = Rmpfr::mpfr(vup/sd, precBits=bits)
    if(!(a<=z &  z<=b)){
        browser()
    }

    ## Separately store and return num&denom of TG
    numer = as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z)))
    denom = as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(a)))

    ## Form p-value as well.
    pv = as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z))/
                    (Rmpfr::pnorm(b)-Rmpfr::pnorm(a)))

    return(list(denom=denom, numer=numer, pv = pv))
}




##' Function to plot qqlot of p-values. Use extra parameter
##' @param pp numeric vector of p-values.
##' @param main label to plot as main title.
qqunif <- function(pp, main=NULL,plot.it=TRUE){
    xy <- stats::qqplot(x=pp,
                 y=seq(from=0,to=1,length=length(pp)), plot.it=FALSE)
    if(plot.it){
        graphics::plot(xy, axes=FALSE)
        graphics::axis(2); graphics::axis(1)
        graphics::abline(0,1)
        if(!is.null(main)) graphics::title(main=main)
    }
    invisible(xy)
}


##' Function to /add/ qq plot points of p-values
##' @param pp numeric vector of p-values.
##' @param main label to plot as main title.
##' @param ... other parameters for \code{qqplot()}.
qqunif_add <- function(pp, main=NULL,...){
    xy <- stats::qqplot(x=pp, y=seq(from=0,to=1,length=length(pp)),plot.it=FALSE)
    graphics::points(xy,...)
}


##' Get all cusums, given start point \code{s} and end point \code{e}
##'
##' @param s Starting index, between \code{1} and \code{n}
##' @param e Ending index, between \code{1} and \code{n}
##' @param y \code{n}-lengthed data vector.
##' @param unsigned \code{TRUE} to return alsolute value.
##' @param simplereturn \code{TRUE} just to return the cusums.
##'
##' @return list of information about the cusum calculations and maximizers in
##'     this interval.
##' @example examples/getcusums-example.R
##' @export
getcusums <- function(s,e,y, unsigned=FALSE, simplereturn=FALSE){

    if(s<=0 | e<= 0) stop("must enter valid e,s >=1 ")

    ## Get all cusum
    cusums =  sapply(s:(e-1), function(b){ cusum(s=s,b=b,e=e,y=y, unsigned=unsigned) })
    if(simplereturn)return(cusums)

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

##' Newer, 10 times faster function for getcusum().
getcusums2 <- function(s, e, cumsums){
    bvec = (s:(e-1))
    n = e-s+1
    cumsums.aug = c(0,cumsums)
    return(-sqrt((e-bvec)/(n*(bvec-s+1)))*(cumsums.aug[bvec+1]-cumsums.aug[s-1+1]) +
        sqrt((bvec-s+1)/(n*(e-bvec)))*(cumsums.aug[e+1]-cumsums.aug[bvec+1]))
}

##' Newer, 10 times faster function for cusum().
cusum2 <- function(s,e,b,cusums){
    cumsums.aug = c(0,cumsums)
    n = e-s+1
    return(-sqrt((e-b)/(n*(b-s+1)))*(cumsums.aug[b+1]-cumsums.aug[s-1+1]) +
        sqrt((b-s+1)/(n*(e-b)))*(cumsums.aug[e+1]-cumsums.aug[b+1]))
}

##' Does a little more than getcusums2().
get_morethan_cusums2 <- function(s,e,cumsums){
    cusums <- getcusums2(s,e,cumsums)
    max.b.ind = which.max(abs(cusums))
    max.b = max.b.ind + s - 1
    max.z = sign(cusums[which.max(cusums)])
    return(list(max.b = max.b,
                max.z = max.z,
                max.cusum = cusums[max.b.ind]))
}



## Either return maximizing breakpoint, or the maximum cusum. If \code{m}
## includes zero, then that is handled to correspond to \code{c(s,e)}
## @param m set of indices that correspond to intervals. If equal to 0, coputes things correspond to \code{c(s,e)}.
## @param s start indices of the interval of interest.
## @param e end index of the interval of interest.
## @param interval set of intervals, produced by \code{generate_intervals()}./
.get_max_b <- function(m,s,e, intervals, y, type=c("cusums","max.b", "max.z")){
    type = match.arg(type)
    if(m!=0){
      s <- intervals$starts[[m]]
      e <- intervals$ends[[m]]
    }
    allcusums = getcusums(s=s,e=e,y=y,unsigned=TRUE)$allcusums
    max.b.before.adding.s <- which.max(allcusums)
    if(type=="cusums"){
        return(max(getcusums(s=s,e=e,y=y,unsigned=TRUE)$allcusums))
    } else if (type == "max.b"){
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


##' Creates bootstrap sample of numeric vector \code{vec}, optionally with a
##' seed.
##' @param vec Numeric vector
##' @param seed Seed for random number generation
##' @return resampled \code{vec}.
##' @export
bootstrap_sample <- function(vec,seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    return(vec[sample.int(length(vec),replace=TRUE)])
}

##' Scale \code{resid} to have empirical std of \code{std}.
##' @param resid numeric vector (meant to be residuals from changepoint model).
##' @param std desired empirical standard deviation
##' @return Properly scaled \code{resid}.
##' @export
scale_resid <- function(resid, std){
    return(resid*(1/stats::sd(resid))*std)
}



##' Only compare things on nodes whose start and end are >1
prune_of_1_length_segments <- function(Tcurr,Scurr,Ecurr){
  Tcurr.copy = Tcurr
  long.enough <- sapply(Tcurr, function(t){
    if(is.null(t)) return(TRUE)
    s = extract(Scurr,t[1],t[2])
    e = extract(Ecurr,t[1],t[2])
    return(e-s>=1)
  })
  Tcurr.copy[!long.enough] = NULL
  return(Tcurr.copy)
}

##' Gets piecewise mean, given segments
piecewise_mean <- function(y,cp){
  stopifnot(all(1<=cp & cp<=length(y)))
  ## stopifnot(all.equal(sort(cp),cp))
  cp = sort(cp)
  segments = lapply(1:(length(cp)+1), function(ii){ v=c(0,cp,length(y));(v[ii]+1):(v[ii+1]) })
  segment.means = sapply(segments, function(mysegment){mean(y[mysegment])})
  cleanmn = rep(NA,length(y))
  lapply(1:length(segments), function(ii){cleanmn[segments[[ii]]] <<- segment.means[ii]})
  return(cleanmn)
}

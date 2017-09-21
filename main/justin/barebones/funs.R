##' Estimation function
##' @param y data vector
estim <- function(y){

    ## code starts here
    n = length(y)
    sums.y = (cumsum(y))
    right.to.left.diff.mn = sapply(1:(n-1), function(i){
        1/(n-i) * (sums.y[n]) - (1/i + 1/(n-i)) * sums.y[i]
    })

    return(which.max(right.to.left.diff.mn))
}


##' polyhedra function
polyhedron <- function(max.cp, n){

    ## Basic checks
    stopifnot(max.cp<n)

    ## Maximizing linear contrast
    max.vec = rep(0,n)
    max.vec[1:max.cp]= -1/max.cp
    max.vec[(max.cp+1):n] = 1/(n-max.cp)

    ## Subtract other contrast from max contrast
    other.cps = (1:(n-1))[-max.cp]
    gamma = t(sapply(other.cps, function(other.cp){
        this.vec = rep(0,n)
        this.vec[1:other.cp] = -1/other.cp
        this.vec[(other.cp+1):n] = +1/(n-other.cp)
        return(max.vec - this.vec)
    }))
    return(gamma)

}


## Make a right-to-left segment contrast
contrast <- function(cp,n){
    vec = rep(0,n)
    vec[1:cp]= -1/cp
    vec[(cp+1):n] = 1/(n-cp)
    return(vec)
}

##' Gets information regarding TG statistic
tg_inf <- function(y, G, u, v, sigma, shift=NULL, nullcontrast=0, bits=50, reduce){

    ## Basic checks
    stopifnot(length(v)==length(y))

    vy = sum(v*y)
    vv = sum(v^2)
    sd = sigma*sqrt(vv)

    ## Shift polyhedron if needed
    if(!is.null(shift)){
        stopifnot(length(shift)==n)
        u = u - G%*%shift
    }


    ## Obtain quantities (intersection points) along the line of v+y
    pvobj <- poly.pval(y=y, G=G, u=u, v=v, sigma=sigma)
    vup = pvobj$vup
    vlo = pvobj$vlo
    vy = max(min(vy, vup),vlo)

    ## Calculate a,b,z for TG = (F(b)-F(z))/(F(b)-F(a))
    z = Rmpfr::mpfr(vy/sd, precBits=bits)
    a = Rmpfr::mpfr(vlo/sd, precBits=bits)
    b = Rmpfr::mpfr(vup/sd, precBits=bits)
    if(!(a<=z &  z<=b)){
        print("F(vlo)<vy<F(vup) was violated, in partition_TG()!")
    }

    ## Separately store and return num&denom of TG
    numer = as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z)))
    denom = as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(a)))

    ## Form p-value as well.
    pv = as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z))/
                    (Rmpfr::pnorm(b)-Rmpfr::pnorm(a)))
    if(!(0 <= pv & pv <= 1)) print("pv was not between 0 and 1, in partition_TG()!")

    return(list(denom=denom, numer=numer, pv = pv, vlo=vlo, vy=vy, vup=vup))
}


##' sampler for RTG
rtg <- function(y, shift, sigma, sigma.add, cp, nsim.inner=100){

    ## Original information
    cp <-  estim(y+shift)
    gamma = polyhedron(cp, n)
    v <- contrast(cp, n)

    ## Get many fudged TG statistics.
    inner.tgs = sapply(1:nsim.inner, function(isim){
        new.noise = rnorm(n,0,sigma.add)
        obj.new = tg_inf(y=y, G=gamma, u=rep(0,nrow(gamma)), shift=new.noise,
                         v=v, sigma=sqrt(sigma^2))
        pv.new = obj.new$pv
        weight.new = obj.new$denom
        return(c(pv.new, weight.new))
    })
    rownames(inner.tgs) = c("pv", "denom")
    pvs = inner.tgs["pv",]
    denoms = inner.tgs["denom",]

    ## Calculate randomized TG statistic
    rtg.pv = sum(pvs*denoms)/sum(denoms)
    ## if(rtg.pv==1) browser()
    return(rtg.pv)
}

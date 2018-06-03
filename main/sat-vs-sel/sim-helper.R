#' Gets selected model p-values.
sel.pvals <- function(y,v,poly,sigma, covariance=NULL, cp.rest, ngen=100){

    ## Generate data
    y0 <- mn + rnorm(n,0,sigma)

    ## Get conditional Gauss parameters
    myparam <- get.cond.gauss.param(cp.rest, y0, covariance=covariance)
    muorig <- myparam$muorig
    Sigmaorig <- myparam$Sigmaorig

    ## Generate new ys, then rejection sample
    enough.samples = FALSE
    accepted=c()
    ys.all = rbind(rep(1,n))[-1,]

    while(!enough.samples){
        ys = (MASS::mvrnorm(n=ngen,mu=muorig,Sigma=Sigmaorig))
        in.polyhedron <- apply(ys, 1, function(myrow){
            return(all(poly$gamma %*% myrow >= poly$u))
        })
        accepted = which(in.polyhedron)
        ys.all = rbind(ys.all,ys[accepted,])
        enough.samples = (nrow(ys.all)>100)
        ## if(!enough.samples) print(nrow(ys.all))
    }

    ## Compute the quantile of the vTY|AY
    vtvec = ys.all %*% cbind(v)

    ## Compute original v^TY
    observed.vt = as.numeric(v%*%y0)

    ## Compute p-value as quantile
    pv = sum(vtvec > observed.vt)/length(vtvec)

    return(pv)
}



##' Get conditional Gaussian parameters, conditioning on the sample means of the
##' conditioned-upon (i.e. most current) model changepoint set.
##' @param cp.curr Current changepoints.
##' @param y0 Observed data vector on hand.
##' @param mn Mean vector, of the observed data vector.
get.cond.gauss.param <- function(cp.curr, y0, mn = rep(0,length(y0)), covariance){


    n=length(y0)

    ## Form the sufficient statistic linear operator
    plateaus = get_plateaus(cp.curr, n)
    suff.rows = do.call(rbind,lapply(plateaus, function(plt){
                              v= rep(0,n); v[plt] = 1/length(plt); return(v)}))
    A = suff.rows

    ## Get null span
    w = A %*% y0
    S = svd(t(A),nu=ncol(A))
    nr = nrow(A)
    A_rest = t(S$u[,(nr+1):ncol(A)])

    Ag = rbind(A,A_rest)

    ## Partition matrix into four blocks
    Sigma = (Ag) %*% covariance %*% t(Ag)
    Si = nrow(A)
    S11 = Sigma[1:Si, 1:Si, drop=FALSE]
    S12 = Sigma[1:Si, (Si+1):n, drop=FALSE]
    S21 = Sigma[(Si+1):n, 1:Si, drop=FALSE]
    S22 = Sigma[(Si+1):n, (Si+1):n, drop=FALSE]

    ## Calculate distr of (A_rest y| Ay) using Schur's complement
    Sigmarest = S22 - S21%*%solve(S11, S12)
    murest = A_rest%*%mn + S21 %*% solve(S11, cbind(w - A %*% mn))

    ## Linear this back to to the original y scale, via left multiplication by Ag
    Aginv = solve(Ag)
    Sigmanew = matrix(0, nrow=n, ncol=n)
    Sigmanew[(Si+1):n, (Si+1):n] = Sigmarest
    Sigmaorig = Aginv %*% Sigmanew %*% t(Aginv)
    muorig = Aginv %*%  cbind(c(w, murest))
    return(list(muorig=muorig, Sigmaorig=Sigmaorig))
}



get_plateaus <- function(inds,n){
    ## Basic checks
    if(length(inds)!=0){
        if(any(is.na(inds))){
            stop("|inds| argument should not contain any NA's! It should be NULL, or equivalently |c()|, if empty.")
        }
    }
    ilist = Map(function(a,b)(a+1):b,
                c(0,inds), c(inds,n))
    return(ilist)
}

##' Do null inference on one-jump data after running fixed number of steps.
dosim <- function(n=200, numSteps=1, lev=0, type=c("gaus", "laplace", "bootstrap"),
                  model = c("binseg", "fusedlasso")){
    mn = c(rep(0, n/2), rep(lev, n/2))
    if(type=="laplace") y = mn + rdoublex(n,lambda=1/sqrt(2))
    if(type=="gaus") y = mn + rnorm(n,0,1)
    if(model=="binseg"){
        h = binSeg_fixedSteps(y, numSteps=numSteps)
    } else {
        h = dualpathSvd2(y=y, D=makeDmat(type="tf",ord=0,m=n), maxsteps=numSteps)
    }
    h.poly = polyhedra(h)
    vlist = make_all_segment_contrasts_from_cp(cp=h$cp, cp.sign=h$cp.sign, n=n)
    retain = which(sapply(vlist, function(v){all.equal(sum(v*mn),0)==TRUE}))
    if(length(retain)==0) return(NULL)
    vlist = vlist[retain]
    pvs = sapply(vlist, function(v){poly.pval2(y=y, v=v, poly=h.poly, sigma=1)$pv})
    return(pvs)
}

##' Do null inference on one-jump data, but stopping with IC.
dosim_icstop <- function(n=200, numSteps=10, lev=0,
                  type=c("gaus", "laplace", "bootstrap"), model = c("binseg", "fusedlasso")){
    mn = c(rep(0, n/2), rep(lev, n/2))
    if(type=="laplace") y = mn + rdoublex(n,lambda=1/sqrt(2))
    if(type=="gaus") y = mn + rnorm(n,0,1)

    if(model=="binseg"){
        h = binSeg_fixedSteps(y, numSteps=numSteps)
    } else {
        h = dualpathSvd2(y=y, D=makeDmat(type="tf",ord=0,m=n), maxsteps=numSteps)
    }
    ic_obj = get_ic(h$cp, h$y, consec=2, sigma=1, type="bic")
    if(ic_obj$flag!="normal")return(NULL)
    vlist = make_all_segment_contrasts_from_cp(cp=h$cp[1:ic_obj$stoptime],
                                               cp.sign=h$cp.sign[1:ic_obj$stoptime],
                                               n=n)
    retain = which(sapply(vlist, function(v){all.equal(sum(v*mn),0)==TRUE}))
    if(length(retain)==0) return(NULL)
    vlist = vlist[retain]
    if(model=="binseg"){
        poly.orig = polyhedra(h, numSteps=ic_obj$stoptime + 2)
    } else {
        h.stopped = dualpathSvd2(y=y, D=makeDmat(type="tf",ord=0,m=n),
                                 maxsteps=ic_obj$stoptime+2)
        poly.orig = polyhedra(h.stopped)
    }
    poly.final = combine.polyhedra(poly.orig, ic_obj$poly)
    pvs = unlist(lapply(vlist, function(v){
        poly.pval2(y=y, v=v, poly=poly.final, sigma=1)$pv
    }))
    return(pvs)
}


dosim_bootstrap <- function(resid, size=length(resid), mn, numSteps=10, 
                            model = c("binseg", "fusedlasso")){
    stopifnot(length(mn)==size)
    resid.boot = bootstrap_sample(vec=resid, size=size)
    y.boot = mn + resid.boot
    if(model=="binseg"){
        h = binSeg_fixedSteps(y.boot, numSteps=numSteps)
    } else {
        h = dualpathSvd2(y=y.boot, D=makeDmat(type="tf",ord=0,m=length(y.boot)), maxsteps=numSteps)
    }
    h.poly = polyhedra(h, numSteps=numSteps)
    vlist = make_all_segment_contrasts_from_cp(cp=h$cp, cp.sign=h$cp.sign, n=length(y.boot))
    retain = which(sapply(vlist, function(v){all.equal(sum(v*mn),0)==TRUE}))
    if(length(retain)==0) return(NULL)
    vlist = vlist[retain]
    pvs = sapply(vlist, function(v){poly.pval2(y=y.boot, v=v, poly=h.poly, sigma=1)$pv})
    return(pvs)
}


dosim_icstop_bootstrap <- function(resid, size=length(resid), mn,
                            model = c("binseg", "fusedlasso")){
    stopifnot(length(mn)==size)
    resid.boot = bootstrap_sample(vec=resid, size=size)
    y.boot = mn + resid.boot
    if(model=="binseg"){
        h = binSeg_fixedSteps(y.boot, numSteps=10)
    } else {
        h = dualpathSvd2(y=y.boot, D=makeDmat(type="tf",ord=0,m=n), maxsteps=10)
    }
    ic_obj = get_ic(h$cp, h$y, consec=2, sigma=1, type="bic")
    if(ic_obj$flag!="normal")return(NULL)
    vlist = make_all_segment_contrasts_from_cp(cp=h$cp[1:ic_obj$stoptime],
                                               cp.sign=h$cp.sign[1:ic_obj$stoptime],
                                               n=length(y.boot))
    retain = which(sapply(vlist, function(v){all.equal(sum(v*mn),0)==TRUE}))
    if(length(retain)==0) return(NULL)
    vlist = vlist[retain]
    if(model=="binseg"){
        poly.orig = polyhedra(h, numSteps=ic_obj$stoptime + 2)
    } else {
        h.stopped = dualpathSvd2(y=y.boot, D=makeDmat(type="tf",ord=0,m=n),
                                 maxsteps=ic_obj$stoptime+2)
        poly.orig = polyhedra(h.stopped)
    }
    poly.final = combine.polyhedra(poly.orig, ic_obj$poly)
    pvs = unlist(lapply(vlist, function(v){
        poly.pval2(y=y.boot, v=v, poly=poly.final, sigma=1)$pv
    }))
    return(pvs)
}

tnorm.surv = function(z, mean=0, sd, a, b) {
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  return((pnorm(b)-pnorm(z))/(pnorm(b)-pnorm(a)))
}

pval = function(y, G, v) {
  z = sum(v*y)
  vv = sum(v^2)
  sd = sqrt(vv)
  
  rho = G %*% v / vv
  vec = (-G %*% y + rho*z) / rho
  quant = G %*% y / G%*%v
  vlo = suppressWarnings(max(vec[rho>0]))
  vlo.ind = which.max(vec[rho>0])
  vup = suppressWarnings(min(vec[rho<0]))
  vty = v%*%y
  ## vty = v%*%y / vv

  vlo.numer = -G %*% y
  vlo.denom = rho * vv

  ## If normal TG, thenn
  p = tnorm.surv(z,0,sd,vlo,vup)
  ## if(plugin) p = 
  return(list(p=p, vlo=vlo, vup=vup, vty=vty, vlo.numer=vlo.numer, vlo.denom=vlo.denom, quant=quant,
              vlo.ind=vlo.ind))
}

## Ryan's function for running 1-step fused lasso vs binary segmentation simulations.
sim = function(n=100, nrep=5000, err.fun=rnorm, mu=rep(0,n), seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  G = matrix(0,2*n-4,n)
  M = A = B = matrix(0,n-1,n)
  for (i in 1:(n-1)) {
    M[i,] = c(rep(-1/i,i),rep(1/(n-i),n-i))
    A[i,] = i*(n-i)/n * M[i,]
    B[i,] = sqrt(i*(n-i)/n) * M[i,]
  }
  
  j.fl = p.fl = vlo.fl = vup.fl = numeric(nrep)
  j.bs = p.bs = vlo.bs = vup.bs = numeric(nrep)
  vty.fl = vty.bs = numeric(nrep)
  maxs.fl = maxs.bs = imaxs.fl = imaxs.bs = list()
  
  for (r in 1:nrep) {
      printprogress(r,nrep)
      y  = mu + err.fun(n)
    
    # Fused lasso
    j1 = which.max(abs(A %*% y))
    s1 = sign(sum(A[j1,]*y))
    G[1:(n-2),] = t(s1*A[j1,] - t(A[-j1,]))
    G[(n-1):(2*n-4),] = t(s1*A[j1,] + t(A[-j1,]))
    obj = pval(y, G, s1*M[j1,])


      ## Added later
      inds = 1:(n-1)
      inds = inds[inds!=j1]
      rownames(G) = c(-inds, inds)
      maxs = order(obj$quant, decreasing=TRUE)[1:10]
      maxs.fl[[r]] = obj$quant[maxs]
      imaxs.fl[[r]] = as.numeric(rownames(G)[maxs])

    j.fl[r] = j1
    p.fl[r] = obj$p
    vlo.fl[r] = obj$vlo
    vup.fl[r] = obj$vup
    vty.fl[r] = obj$vty

    # Binary segmentation
    j2 = which.max(abs(B %*% y))
    s2 = sign(sum(B[j2,]*y))
    G[1:(n-2),] = t(s2*B[j2,] - t(B[-j2,]))
    G[(n-1):(2*n-4),] = t(s2*B[j2,] + t(B[-j2,]))
    obj = pval(y, G, s2*M[j2,])

      ## Added later
      inds = 1:(n-1)
      inds = inds[inds!=j2]
      rownames(G) = c(-inds, inds)
      maxs = order(obj$quant, decreasing=TRUE)[1:10]
      maxs.bs[[r]] = obj$quant[maxs]
      imaxs.bs[[r]] = as.numeric(rownames(G)[maxs])

    j.bs[r] = j2
    p.bs[r] = obj$p
    vlo.bs[r] = obj$vlo
    vup.bs[r] = obj$vup
    vty.bs[r] = obj$vty
  }
  
  # Changepoint histograms
  ## pdf(file="~/Desktop/hist.pdf",width=10,height=5)
  ylim = range(c(hist(j.fl, breaks=0:(n-1), plot=FALSE)$density,
                 hist(j.bs, breaks=0:(n-1), plot=FALSE)$density))
  hist(j.fl, col="pink", breaks=0:(n-1), prob=TRUE, ylim=ylim,
       xlab="Jump location", main="Histogram of jump locations")
  hist(j.bs, col=rgb(0,0.5,0.5,0.5), breaks=0:(n-1), prob=TRUE, add=TRUE)
  legend("topright", legend=c("FL","BS"), col=c("pink",rgb(0,0.5,0.5)), lty=1)
  ## graphics.off()

  # QQ plots
  ## pdf(file="~/Desktop/qq.pdf",width=10,height=5)
  plot(1:nrep/nrep, 1:nrep/nrep, type="l", lty=2, col="gray", 
       xlab="Expected", ylab="Observed")
  points(1:nrep/nrep, sort(p.fl), col="pink", pch=20)
  points(1:nrep/nrep, sort(p.bs), col=rgb(0,0.5,0.5,0.5), pch=20)
  legend("topright", legend=c("FL","BS"), col=c("pink",rgb(0,0.5,0.5)), pch=20)
  ## graphics.off()
  
  invisible(list(j.fl=j.fl, p.fl=p.fl, vlo.fl=vlo.fl, vup.fl=vup.fl,
                 j.bs=j.bs, p.bs=p.bs, vlo.bs=vlo.bs, vup.bs=vup.bs,
                 vty.bs=vty.bs, vty.fl=vty.fl, maxs.bs=maxs.bs, imaxs.bs=imaxs.bs, maxs.fl=maxs.fl,
                 imaxs.fl=imaxs.fl))
}



## Helpers
get_sbe <- function(myrow){
    ## print(myrow)
    tol = 1E-10
    myrow[abs(myrow)<tol] = 0
    v = (makeDmat(length(myrow),order=0) %*% myrow)
    v[abs(v)<tol] = 0
    s = min(which(myrow!=0))
    e = max(which(myrow!=0))
    myrow.sub = myrow[s:e]
    vsub = (makeDmat(length(myrow.sub),order=0) %*% myrow.sub)
    vsub[abs(vsub)<tol] = 0
    b = which(vsub!=0) + s - 1
    if(length(b)>1) browser()
    list(s=s,
         b=b,
         e=e)
}
make.base.contrast.from.sbe <- function(s,b,e,n){
    L = s:b
    R = (b+1):e
    v = rep(0,n)
    ## v[L] = 1/length(L)
    v[L] = -1/length(L)
    v[R] = 1/length(R)
    return(v)
}

## item-by-item comparison of the things that are maximized in forming Vlo. |y| is assumed to be given
vlo.surgery = function(y) {

    n = length(y)
  
  ## Convenient quantities, formed once
  G = matrix(0,2*n-4,n)
  M = A = B = matrix(0,n-1,n)
  for (i in 1:(n-1)) {
    M[i,] = c(rep(-1/i,i),rep(1/(n-i),n-i))
    A[i,] = i*(n-i)/n * M[i,]
    B[i,] = sqrt(i*(n-i)/n) * M[i,]
  }
  
  ## vlo.denom.fl = vlo.numer.fl = vlo.denom.bs = vlo.numer.bs = list()

    ## y  = mu + err.fun(n)
    
    # Fused lasso
    j1 = which.max(abs(A %*% y))
    s1 = sign(sum(A[j1,]*y))
    G[1:(n-2),] = t(s1*A[j1,] - t(A[-j1,]))
    G[(n-1):(2*n-4),] = t(s1*A[j1,] + t(A[-j1,]))
    obj = pval(y, G, s1*M[j1,])
    j.fl = j1
    p.fl = obj$p

    vlo.fl = obj$vlo
    vlo.denom.fl = obj$vlo.denom
    vlo.numer.fl = obj$vlo.numer

    vup.fl = obj$vup
    vty.fl = obj$vty

    # Binary segmentation
    j2 = which.max(abs(B %*% y))
    s2 = sign(sum(B[j2,]*y))
    G[1:(n-2),] = t(s2*B[j2,] - t(B[-j2,]))
    G[(n-1):(2*n-4),] = t(s2*B[j2,] + t(B[-j2,]))
    obj = pval(y, G, s2*M[j2,])
    j.bs = j2
    p.bs = obj$p
    vlo.bs = obj$vlo

    vlo.denom.bs = obj$vlo.denom
    vlo.numer.bs = obj$vlo.numer

    vup.bs = obj$vup
    vty.bs = obj$vty

  invisible(list(j.fl=j.fl,
                 p.fl=p.fl,
                 vlo.fl=vlo.fl,
                 vup.fl=vup.fl,
                 j.bs=j.bs,
                 p.bs=p.bs,
                 vlo.bs=vlo.bs,
                 vup.bs=vup.bs,
                 vty.bs=vty.bs,
                 vty.fl=vty.fl,
                 vlo.denom.fl = vlo.denom.fl,
                 vlo.numer.fl =  vlo.numer.fl,
                 vlo.denom.bs = vlo.denom.bs,
                 vlo.numer.bs  = vlo.numer.bs,
                 v=v,
                 ))
}




## Modifying Ryan's function for running 1-step fused lasso vs binary
## segmentation simulations.
simboot <- function(n=100, nrep=5000, err.fun=rnorm, mu=rep(0,n), seed=NULL, sigma=1, samp=NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  G = matrix(0,2*n-4,n)
  M = A = B = matrix(0,n-1,n)
  for (i in 1:(n-1)) {
    M[i,] = c(rep(-1/i,i),rep(1/(n-i),n-i))
    A[i,] = i*(n-i)/n * M[i,]
    B[i,] = sqrt(i*(n-i)/n) * M[i,]
  }
  
  j.fl = p.fl = vlo.fl = vup.fl = numeric(nrep)
  j.bs = p.bs = vlo.bs = vup.bs = numeric(nrep)
  vty.fl = vty.bs = numeric(nrep)
  p2.bs = p2.fl = numeric(nrep)

  
  start.time=Sys.time()
  for (r in 1:nrep) {
      printprogress(r,nrep, start.time=start.time)

      resids = err.fun(n, samp=samp)
      y  = mu + resids

      ## Added: matrix of bootstrapped residuals.
      nboot = 3000
      residmat = t(sapply(1:nboot, function(iboot){
          resids[sample(n, size=n, replace=TRUE)]
      }))
    
    # Fused lasso
    j1 = which.max(abs(A %*% y))
    s1 = sign(sum(A[j1,]*y))
    G[1:(n-2),] = t(s1*A[j1,] - t(A[-j1,]))
    G[(n-1):(2*n-4),] = t(s1*A[j1,] + t(A[-j1,]))
    obj = binSegInf::poly.pval(y=y, G=G, v=s1*M[j1,], u=rep(0,nrow(G)), sigma=sigma)
    p2 = pval_plugin(obj$vlo, obj$vup, obj$vty, s1*M[j1,], bootmat=residmat) ##added
    j.fl[r] = j1
    p.fl[r] = obj$p
    p2.fl[r] = p2 ## Added
    vlo.fl[r] = obj$vlo
    vup.fl[r] = obj$vup
    vty.fl[r] = obj$vty

    # Binary segmentation
    j2 = which.max(abs(B %*% y))
    s2 = sign(sum(B[j2,]*y))
    G[1:(n-2),] = t(s2*B[j2,] - t(B[-j2,]))
    G[(n-1):(2*n-4),] = t(s2*B[j2,] + t(B[-j2,]))
    obj = binSegInf::poly.pval(y=y, G=G, v=s2*M[j2,], u=rep(0,nrow(G)), sigma=sigma)
    p2 = pval_plugin(obj$vlo, obj$vup, obj$vty, s2*M[j2,], bootmat=residmat) ##added
    j.bs[r] = j2
    p.bs[r] = obj$p
    p2.bs[r] = p2 ## Added
    vlo.bs[r] = obj$vlo
    vup.bs[r] = obj$vup
    vty.bs[r] = obj$vty
  }
  
  # Changepoint histograms
  ## pdf(file="~/Desktop/hist.pdf",width=10,height=5)
  ## ylim = range(c(hist(j.fl, breaks=0:(n-1), plot=FALSE)$density,
  ##                hist(j.bs, breaks=0:(n-1), plot=FALSE)$density))
  ## hist(j.fl, col="pink", breaks=0:(n-1), prob=TRUE, ylim=ylim,
  ##      xlab="Jump location", main="Histogram of jump locations")
  ## hist(j.bs, col=rgb(0,0.5,0.5,0.5), breaks=0:(n-1), prob=TRUE, add=TRUE)
  ## legend("topright", legend=c("FL","BS"), col=c("pink",rgb(0,0.5,0.5)), lty=1)
  ## graphics.off()

  # QQ plots
  ## pdf(file="~/Desktop/qq.pdf",width=10,height=5)
  ## plot(1:nrep/nrep, 1:nrep/nrep, type="l", lty=2, col="gray", 
  ##      xlab="Expected", ylab="Observed")
  ## points(1:nrep/nrep, sort(p.fl), col="pink", pch=20)
  ## points(1:nrep/nrep, sort(p.bs), col=rgb(0,0.5,0.5,0.5), pch=20)
  ## legend("topright", legend=c("FL","BS"), col=c("pink",rgb(0,0.5,0.5)), pch=20)
  ## ## graphics.off()
  
  invisible(list(j.fl=j.fl, p.fl=p.fl, vlo.fl=vlo.fl, vup.fl=vup.fl,
                 j.bs=j.bs, p.bs=p.bs, vlo.bs=vlo.bs, vup.bs=vup.bs,
                 vty.bs=vty.bs, vty.fl=vty.fl,
                 p.bs=p.bs,
                 p.fl=p.fl,
                 p.bootstrap.bs=p2.bs,
                 p.bootstrap.fl=p2.fl
                 ))
}

simnew = simboot



## Helpers
get_sbe <- function(myrow){
    ## print(myrow)
    tol = 1E-10
    myrow[abs(myrow)<tol] = 0
    v = (makeDmat(length(myrow),order=0) %*% myrow)
    v[abs(v)<tol] = 0
    s = min(which(myrow!=0))
    e = max(which(myrow!=0))
    myrow.sub = myrow[s:e]
    vsub = (makeDmat(length(myrow.sub),order=0) %*% myrow.sub)
    vsub[abs(vsub)<tol] = 0
    b = which(vsub!=0) + s - 1
    if(length(b)>1) browser()
    list(s=s,
         b=b,
         e=e)
}
make.base.contrast.from.sbe <- function(s,b,e,n){
    L = s:b
    R = (b+1):e
    v = rep(0,n)
    ## v[L] = 1/length(L)
    v[L] = -1/length(L)
    v[R] = 1/length(R)
    return(v)
}

## item-by-item comparison of the things that are maximized in forming Vlo. |y| is assumed to be given
vlo.surgery = function(y) {

    n = length(y)
  
  ## Convenient quantities, formed once
  G = matrix(0,2*n-4,n)
  M = A = B = matrix(0,n-1,n)
  for (i in 1:(n-1)) {
    M[i,] = c(rep(-1/i,i),rep(1/(n-i),n-i))
    A[i,] = i*(n-i)/n * M[i,]
    B[i,] = sqrt(i*(n-i)/n) * M[i,]
  }
  
  ## vlo.denom.fl = vlo.numer.fl = vlo.denom.bs = vlo.numer.bs = list()

    ## y  = mu + err.fun(n)
    
    # Fused lasso
    j1 = which.max(abs(A %*% y))
    s1 = sign(sum(A[j1,]*y))
    G[1:(n-2),] = t(s1*A[j1,] - t(A[-j1,]))
    G[(n-1):(2*n-4),] = t(s1*A[j1,] + t(A[-j1,]))
    obj = pval(y, G, s1*M[j1,])
    j.fl = j1
    p.fl = obj$p

    vlo.fl = obj$vlo
    vlo.denom.fl = obj$vlo.denom
    vlo.numer.fl = obj$vlo.numer

    vup.fl = obj$vup
    vty.fl = obj$vty

 ## Binary segmentation
    j2 = which.max(abs(B %*% y))
    s2 = sign(sum(B[j2,]*y))
    G[1:(n-2),] = t(s2*B[j2,] - t(B[-j2,]))
    G[(n-1):(2*n-4),] = t(s2*B[j2,] + t(B[-j2,]))
    obj = pval(y, G, s2*M[j2,])
    j.bs = j2
    p.bs = obj$p
    vlo.bs = obj$vlo

    vlo.denom.bs = obj$vlo.denom
    vlo.numer.bs = obj$vlo.numer

    vup.bs = obj$vup
    vty.bs = obj$vty

  invisible(list(j.fl=j.fl,
                 p.fl=p.fl,
                 vlo.fl=vlo.fl,
                 vup.fl=vup.fl,
                 j.bs=j.bs,
                 p.bs=p.bs,
                 vlo.bs=vlo.bs,
                 vup.bs=vup.bs,
                 vty.bs=vty.bs,
                 vty.fl=vty.fl,
                 vlo.denom.fl = vlo.denom.fl,
                 vlo.numer.fl =  vlo.numer.fl,
                 vlo.denom.bs = vlo.denom.bs,
                 vlo.numer.bs  = vlo.numer.bs,
                 v=v,
                 ))
}


## A set of error functions to use for simulations.
lapl <- function(n,samp=NULL){ rexp(n,rate=sqrt(2)) * sample(c(-1,1),n,replace=TRUE)}
rt2 <- function(n,samp=NULL){ rt(n, df=2) }
rt3 <- function(n,samp=NULL,scale=FALSE){
    df = 3
    if(scale){ std = sqrt(df/(df-2)) } else { std = 1 }
    return(rt(n, df=3)/std)
}
## rresid <- function(resid, sigma){
##     n = length(resid)
##     bootstrap.inds = bootstrap_ind(n, size=n)
##     return(resid[bootstrap.inds])
## }

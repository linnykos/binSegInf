##' Run simulations for a single setting for SBS, for the one-jump
##' signal.
do.one.sim = function(delta,nsim,sigma,numsteps,n=12){
  p = rep(NA,nsim)
  theta = c(rep(0,n/2), rep(delta,n/2))
  for(isim in 1:nsim){
    y = theta + rnorm(n,0,sigma)
    a = binseg.by.size(y, numsteps,verbose=FALSE)
    if(a$B[1]==n/2){
      v = make.v(a$B[1],a$B,a$Z,n)
      p[isim]=poly.pval(y,a$G,a$u,v,sigma)$pv
    }
  }
  return(p)
}


##' Run simulations for a single setting for 1d fused lasso, for the
##' one-jump signal.
do.one.sim.1dfl = function(delta,nsim,sigma,numsteps,n=12){
  p = rep(NA,nsim)
  theta = c(rep(0,n/2), rep(delta,n/2))
  for(isim in 1:nsim){
    y = theta + rnorm(n,0,sigma)
    path   = dualpathSvd2(y, dual1d_Dmat(n), maxsteps = numsteps, approx = TRUE)
    G      = path$Gammat[1:path$nk[1],]
    d      = getdvec(obj=path, y=y0, k=1, usage = "dualpathSvd", type=testtype, matchstep=TRUE)
    if(path$pathobj$B[1] == n/2){
        pvals.correct[jj] = pval.fl1d(y=y0, G=G, dik=d, sigma=sigma, approx=TRUE, threshold=TRUE, approxtype="rob")
        jj = jj + 1
        if(verbose) cat("\r", jj, "of", nsim)
      }


    if(a$B[1]==n/2){
      v = make.v(a$B[1],a$B,a$Z,n)
      p[isim]=poly.pval(y,a$G,a$u,v,sigma)$pv
    }
  }
  return(p)
}

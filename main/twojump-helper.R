
## Define some helper functions
do.one.sim.twojump = function(delta, nsim, sigma, numsteps, nstepmodel, n=12){
  
  stopifnot(numsteps==nstepmodel)
  
  produce.pval.onestep = function(a,y){
    if(a$B[1] ==n/3 ){
      v = make.v(a$B[1],a$B,a$Z,n)
      return(poly.pval(y,a$G,a$u,v,sigma)$pv)
    } else {
      return(NA)
    }
  }
  
  produce.pval.twostep = function(a,y){
    if(a$B[1]==n/3 && a$B[2] == 2*n/3){
      v = make.v(a$B[1],a$B,a$Z,n)
      return(poly.pval(y,a$G,a$u,v,sigma)$pv)
    } else {
      return(NA)
    }
  }
  
  produce.pval = (if(nstepmodel==1) produce.pval.onestep else produce.pval.twostep)
  
  
  ## Simulation settings
  p = rep(NA,nsim)
  theta = c(rep(0,n/3), 
            rep(delta,n/3),
            rep(0,n/3))
  
  ## Run simulations
  for(isim in 1:nsim){
    y = theta + rnorm(n,0,sigma)
    a = binseg.by.size(y, numsteps,verbose=FALSE)
    p[isim] = produce.pval(a,y)
  }
  return(p)
}



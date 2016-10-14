do.one.sim.threejump = function(delta, nsim, sigma, numsteps=1, n=16){
  ## Simulation settings
  theta = c(rep(0,n/4), 
            rep(delta,n/4),
            rep(-delta,n/4),
            rep(0,n/4))
  correct.locs = c(n/4, n/2, 3*n/4)
  
    
  ## Run simulations
  locs = rep(NA,nsim)
  for(isim in 1:nsim){
    y = theta + rnorm(n,0,sigma)
    a = binseg.by.size(y, numsteps,verbose=FALSE)
    locs[isim] = a$B[numsteps]
  }
  return(locs)
}
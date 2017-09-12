## Synopsis: We want to see uniform null p-values for fused lasso
## Simulation settings.
source("../main/justin/sim-helper.R")
source("../main/justin/sim-driver.R")
numSteps = 1
nsim = 100
nsim.is = 100
lev = 0
## pvec = rep(NA,nsim)

## cat(fill=TRUE)
## pvec = mclapply(1:nsim, function(isim){
##     cat('\r', isim, "out of", nsim)
##     my.mn <- meanfun(lev,n)
##     y = (my.mn + stats::rnorm(n,0,sigma))

##     ## Do fused lasso inference
##     D = makeDmat(n,type='tf',ord=0)
##     obj <- genlassoinf::dualpathSvd2(y,D=D,maxsteps=numSteps,approx=TRUE)
##     contrasts <- make_all_segment_contrasts(obj)
##     pvec = pvec.plain = setNames(rep(NA,length(obj$cp)), obj$cp)
##     poly <- polyhedra(obj$Gobj.naive$G, obj$Gobj.naive$u)

##     ## pv <- poly.pval2(y=y, poly=poly, v=contrasts[[ii]], sigma=sigma,
##     ##                              reduce=reduce)$pv
##     pv = randomized_genlasso_pv(y=y, v=contrasts[[ii]],
##                                        sigma=sigma, numSteps=numSteps,
##                                        numIntervals=numIntervals,
##                                        nsim.is=nsim.is, bits=100,
##                                        reduce=reduce, augment=augment)
##     return(pv)

## }, mc.cores = 3)
## hist(unlist(pvec),xlim=c(0,1))



    ## sigma=1; lev=0; nsim.is=100; numSteps=1;
    ## numIntervals=20; n=6; meanfun=onejump;
    ## reduce=TRUE;augment=TRUE;  bootstrap=FALSE; std.bootstrap=NULL;
    ## cleanmn.bootstrap=NULL;
    ## type = "randomized"

## What about WBS?
sim.settings = list(sigma=1, lev=0, nsim.is=100, numSteps=1,
                    numIntervals=10, n=4, meanfun=onejump,
                    reduce= FALSE, augment=TRUE, bootstrap=FALSE, std.bootstrap=NULL,
                    cleanmn.bootstrap=NULL,
                    type = "randomized")

## pv = onesim_wbs(sim.settings)

nsim = 100
cat(fill=TRUE)
pvec = unlist(mclapply(1:nsim, function(isim){cat('\r',isim,'out of',nsim);onesim_wbs(sim.settings)}, mc.cores=3))
## pvec.plain = unlist(mclapply(1:nsim, function(isim){cat('\r',isim,'out of',nsim);onesim_wbs(sim.settings)}, mc.cores=3))
## pvec.naive = unlist(mclapply(1:nsim, function(isim){cat('\r',isim,'out of',nsim);onesim_naive(sim.settings)}, mc.cores=3))
## source("../main/justin/sim-driver.R")
## a = onesim_naive(sim.settings)

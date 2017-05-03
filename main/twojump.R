## Simulation settings
n = 30
nsim.is = 100
numSteps = 2
numIntervals = 10
n.levs = 4
levs=seq(from=0,to=3,length=n.levs)
nsims = seq(from=1000,to=500,length=n.levs)
sigma = 1
mc.cores = 8
sim.settings <- list(levs = levs,
                     nsim.is = nsim.is,
                     numSteps = numSteps,
                     numIntervals = numIntervals,
                     n = n,
                     mn = mn.twojump,
                     nsims = nsims,
                     sigma = sigma,
                     bootstrap = FALSE)

sim_driver(sim.settings=sim.settings,
           filename="twojump.Rdata",
           dir="../results",
           mc.cores=mc.cores)




## ## First check
## pmat.wbsfs.plain = matrix(NA,nrow=500, ncol=15)
## for(isim in 1:500){
##     print(isim)

##     my.mn <- mn(lev,n)
##     set.seed(isim)
##     y <- my.mn + rnorm(n,0,sigma)
##     ## Why are p-values weird
##     method <- wildBinSeg_fixedSteps
##     set.seed(isim+500)
##     intervals <- generate_intervals(n=length(y),numIntervals=numIntervals)
##     obj <- method(y, numSteps=numSteps, intervals=intervals)
##     contrast <- make_all_segment_contrasts(obj)
##     poly <- polyhedra(obj, v=v, reduce=TRUE)
##     for(ii in 1:length(obj$cp)){
##         pmat.wbsfs.plain[isim,obj$cp[ii]] <- poly.pval(y=y, G=poly$gamma, u=poly$u,
##                                                     v=contrast[[ii]], sigma=sigma)$pv
##     }
## }
## qqunif(pmat.wbsfs.plain[,5],ylim=c(0,1),xlim=c(0,1))
## qqunif_add(pmat.wbsfs.plain[,10],col='red')

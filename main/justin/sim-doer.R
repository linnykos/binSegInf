## Synopsis: Actually run some simulations
source("../main/justin/sim-driver.R")
sim.settings = list(sigma=1, lev=0, nsim.is=50, numSteps=1,
                    numIntervals=5, n=10, meanfun=onejump,
                    reduce=FALSE,augment=TRUE,  bootstrap=FALSE, std.bootstrap=NULL,
                    cleanmn.bootstrap=NULL, thresh = 1,
                    type = "random",v = runif(10,-1,1))##plain
sim.settings.plain = sim.settings; sim.settings.plain[["type"]]="plain"
nsim=10
onesim_wbs(sim.settings)
ps = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_wbs(sim.settings)}, mc.cores=3)


## Trying the winning-est interval idea.
source("../main/justin/sim-driver.R")
nsim=1000
pvs = mclapply(1:nsim, function(isim){printprogress(isim, nsim);onesim_wbs_new()}, mc.cores=3)
qqunif(unlist(pvs))


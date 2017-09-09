## Simulation settings
n = 30
nsim.is = 100
numSteps = 2
numIntervals = 10
n.levs = 2
levs=seq(from=0,to=2,length=n.levs)
nsims = seq(from=1000,to=500,length=n.levs)/2
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
                     ## resid.cleanmn = resid.cleanmn)



sim_driver(sim.settings=sim.settings,
           filename="twojump.Rdata",
           dir="../results",
           mc.cores=mc.cores)



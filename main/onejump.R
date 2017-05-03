## Simulation settings
n = 20
nsim.is = 100
numSteps = 1
numIntervals = 20
n.levs = 5
levs = seq(from=0,to=2,length=n.levs)
nsims = seq(from=1000,to=500,length=n.levs)
sigma = 1
sim.settings <- list(levs = levs,
                     nsim.is = nsim.is,
                     numSteps = numSteps,
                     numIntervals = numIntervals,
                     n = n,
                     mn = mn.onejump,
                     nsims = nsims,
                     sigma = sigma)

sim_driver(sim.settings,"onejump.Rdata", dir="../results", reduce=TRUE)


sim_driver(sim.settings=sim.settings,
           filename="onejump.Rdata",
           dir="../results",
           mc.cores=mc.cores)

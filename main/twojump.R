## Simulation settings
n = 30
nsim.is = 100
numSteps = 2
numIntervals = 10
n.levs = 2
levs=seq(from=0,to=2,length=n.levs)
nsims = seq(from=1000,to=500,length=n.levs)
sigma = 1
sim.settings <- list(levs = levs,
                     nsim.is = nsim.is,
                     numSteps = numSteps,
                     numIntervals = numIntervals,
                     n = n,
                     mn = mn.twojump,
                     nsims = nsims,
                     sigma = sigma)

sim_driver(sim.settings,"twojump.Rdata", dir="../data")

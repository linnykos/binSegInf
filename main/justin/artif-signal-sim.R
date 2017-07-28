## Load in data and simulation driver functions
load("../data/coriell.Rdata")
source('../main/justin/sim-driver.R')

## n = length(coriell_mn(1))
n = 10
numSteps = 1
numIntervals = 100
lev=3
sigma=1
nreplicate = 100
mc.cores=3
sim.settings = list(numIntervals=numIntervals,
                    nreplicate=nreplicate,
                    lev=lev,
                    numSteps=numSteps,
                    nsim.is=nsim.is)
nsim.is = 100

## Generate data
set.seed(0)
## y <- coriell_mn(lev,n) + rnorm(n,0,std)
y <- onejump(lev,n) + rnorm(n,0,1)

nreplicate = 10


## Simulation settings
sim.settings = list(sigma=1, lev=1, nsim.is=100, numSteps=1,
                    numIntervals=100, n=10, meanfun=onejump,
                    reduce=FALSE,augment=TRUE,  bootstrap=FALSE, std.bootstrap=NULL,
                    cleanmn.bootstrap=NULL,
                    type = "random")##plain

## Actually run the simulations
## a = onesim_sbs(sim.settings)
a = onesim_wbs(sim.settings)
a = onesim_fusedlasso(sim.settings)

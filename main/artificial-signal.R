load("../data/coriell.Rdata")

##' Takes a given mean and multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1,n){
    newmn2 = (newmn[1101:1300][seq(from=1,to=200,length=100)])
    h = max(abs(newmn))
    return((newmn2 / h * std) * lev)
}
n = length(coriell_mn(1))
nsim.is = 10
numSteps = 5
numIntervals = 30
n.levs = 2
levs=seq(from=0,to=3,length=n.levs)
## nsims = seq(from=100,to=50,length=n.levs)#rep(1,n.levs)#
nsims =rep(1,n.levs) 
sim.settings <- list(levs = levs,
                     nsim.is = nsim.is,
                     numSteps = numSteps,
                     numIntervals = numIntervals,
                     n = n,
                     mn = coriell_mn,
                     nsims = nsims,
                     sigma = std)

sim_driver(sim.settings,"artificial-smaller.Rdata", dir="../data")

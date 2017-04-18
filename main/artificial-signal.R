load("data/coriell.Rdata")

##' Takes a given mean and multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1,n){
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}
n = length(coriell_mn(1))
nsim.is = 1
numSteps = 5
numIntervals = 10
n.levs = 1
levs=seq(from=3,to=3,length=n.levs)
nsims = rep(1,n.levs)#seq(from=1000,to=500,length=n.levs)
sim.settings <- list(levs = levs,
                     nsim.is = nsim.is,
                     numSteps = numSteps,
                     numIntervals = numIntervals,
                     n = n,
                     mn = coriell_mn,
                     nsims = nsims,
                     sigma = std)


sim_driver(sim.settings,"artificial.Rdata", dir="../data")

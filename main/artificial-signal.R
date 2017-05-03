load("../data/coriell.Rdata")

##' Takes a given mean and multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1,n){
    ## newmn = (newmn[1101:1300][seq(from=1,to=200,length=100)])
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}
n = length(coriell_mn(1))
nsim.is = 1
numSteps = 5
numIntervals = 50
n.levs = 1
levs=seq(from=0,to=3,length=n.levs)
## nsims = seq(from=100,to=50,length=n.levs)#rep(1,n.levs)#
nsims = 1
nsims =rep(1,n.levs)

## Seeing the memory
method <- wildBinSeg_fixedSteps
intervals <- generate_intervals(n=length(y),numIntervals=numIntervals)
obj <- method(y, numSteps=numSteps, intervals=intervals)
contrast <- make_all_segment_contrasts(obj)
poly <- polyhedra(obj,v=v,reduce=reduce)



sim.settings <- list(levs = levs,
                     nsim.is = nsim.is,
                     numSteps = numSteps,
                     numIntervals = numIntervals,
                     n = n,
                     mn = coriell_mn,
                     nsims = nsims,
                     sigma = std)

sim_driver(sim.settings = sim.settings,
           filename = "artificial.Rdata",
           dir="../data",
           reduce=TRUE)

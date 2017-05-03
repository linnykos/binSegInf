load("../data/coriell.Rdata")

##' Takes a given mean and multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1,n){
    ## newmn = (newmn[1101:1300][seq(from=1,to=200,length=100)])
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}

n = length(coriell_mn(1))
nsim.is = 500
numSteps = 5
numIntervals = 500
n.levs = 1
levs=seq(from=0,to=3,length=n.levs)
nsims = seq(from=100,to=50,length=n.levs)
bootstrap=TRUE
reduce=TRUE

## sim.settings <- list(levs = levs,
##                      nsim.is = nsim.is,
##                      numSteps = numSteps,
##                      numIntervals = numIntervals,
##                      n = n,
##                      mn = coriell_mn,
##                      nsims = nsims,
##                      sigma = std,
##                      std = std,
##                      bootstrap=bootstrap
##                      )


## sim_driver(sim.settings = sim.settings,
##            ## filenames = paste0("artificial-lev-",mylev,".Rdata"),
##            filenames = "artificial.Rdata",
##            dir="../data", reduce=reduce )


lev=4
mn <- coriell_mn
set.seed(0)
y <- mn(lev,n) + rnorm(n,0,std)
method <- wildBinSeg_fixedSteps
intervals <- generate_intervals(n=length(y),numIntervals=numIntervals)

print("This is how long it took to run the algorithm:")
mytime <- proc.time()
obj <- method(y, numSteps=numSteps, intervals=intervals, verbose=TRUE)
mynewtime <- proc.time()
print(mynewtime - mytime)

contrast <- make_all_segment_contrasts(obj)
v = contrast[[1]]

print("This is how long it took to make the polyhedra:")
mytime <- proc.time()
poly <- polyhedra(obj,v=v,reduce=reduce)
mynewtime <- proc.time()
print(mynewtime - mytime)
print(format(object.size(poly), "Mb"))

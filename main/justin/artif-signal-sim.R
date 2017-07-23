load("../data/coriell.Rdata")

##' Takes a given mean and multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1,n){
    ## newmn2 = (newmn[1101:1300][seq(from=1,to=200,length=100)])
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}


mn.onejump <- function(lev=1,n){
    c(rep(0,n/2),rep(lev,n/2))
}
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
y <- mn.onejump(lev,n) + rnorm(n,0,1)

p.wbsfs.list = mclapply(1:nreplicate, function(irep){
    cat("\r", "replicate", irep, "out of", nreplicate[ii])

    method = wildBinSeg_fixedSteps
    intervals = generate_intervals(length(y), numIntervals)
    obj = method(y, numSteps=numSteps, intervals=intervals)
    poly = polyhedra(obj)

    vs = make_all_segment_contrasts(obj)
    p.wbsfs = rep(NA,length(obj$cp))
    names(p.wbsfs) = obj$cp
    for(ii in 1:length(obj$cp)){
        p.wbsfs[ii] <- randomized_wildBinSeg_pv(y=y, v=vs[[ii]], sigma=sigma,
                                                numSteps=numSteps,
                                                numIntervals=numIntervals,
                                                nsim.is=nsim.is, bits=100,
                                                augment=TRUE)
    }
    return(p.wbsfs)
}, mc.cores=mc.cores)

## Save result
filename = paste0("artif.Rdata")
save(p.wbsfs.list, sim.settings, file = file.path("../../results", filename))

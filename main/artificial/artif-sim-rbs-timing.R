# Synopsis: See how long simulations are going to take.
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))
library(microbenchmark)

##' Multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1, n){
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}
fac=2
sigma=1
sigma = sd(y.orig[1:200]) * fac
sigma.add = sigma*0.2
set.seed(0)
y = newmn[-(1:200)] + bootstrap_sample(resid.cleanmn[-(1:200)]) * fac

## Timing the inference.
verbose=TRUE
bits=1000
numIS=200
min.num.things=100
how.close=5
mytime = list()
pvs.rbs = c()
for(ii in 1:5){
    mytime[[ii]] = microbenchmark({
        start.time = Sys.time()
        pvs.rbs[[ii]] = do_rbs_inference(y=y, max.numSteps=10, consec=2, sigma=sigma,
                                       postprocess=TRUE, locs=1:length(y), numIS=numIS, min.num.things=min.num.things,
                                       inference.type="pre-multiply", bits=bits, sigma.add=sigma.add, verbose=verbose,
                                       start.time=start.time,
                                       whichv=ii,
                                       how.close=how.close, return.more.things=TRUE)
        cat(fill=TRUE)
    },times=1)
}

length(pvs.rbs)
results = pvs.rbs
## First factor
for(ii in 1:3){
    results = results.by.fac[[ii]][[1]]
    ## First p-value
    pvs = result$parts.so.far
    weights = result$parts.so.far
}


## No effect of time; 10 minutes each.

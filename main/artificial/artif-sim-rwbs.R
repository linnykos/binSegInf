# Data directory
datadir = "../data"
outputdir = "../output"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))

## Simulation driver
onesim_rwbs <- function(y.orig){

    ## Add bootstrapped residuals around a cleaned mean, with known sigma
    ## estimated in-sample.
    sigma = sd(y.orig[1:200])
    resid = resid.cleanmn[-(1:200)]
    resid.order = bootstrap_ind(length(resid))
    y = newmn[-(1:200)] + resid[resid.order]

    ## The idea is to add bootstrap noise, then see the (1) conditional and
    ## unconditional powers of the tests done at the viscinities of the true
    ## guys, and (2) uniformity (or lack thereof) in the tests conducted
    ## regarding null locations.

    lev=3
    y = rnorm(10,0,1) + c(rep(0,5), rep(lev,5))
    sigma=1

    ## Just do two types of inference: randwbs, randBS
    intervals = intervals(numIntervals=length(y), n=length(y))
source(file=file.path("../main/artificial/artif-helpers.R"))
    results.rwbs = do_rwbs_inference(y=y, max.numSteps=10, numIntervals=NULL,
                                 intervals=intervals, consec=2, sigma=sigma,
                                 postprocess=TRUE, better.segment=TRUE,
                                 locs=1:length(y), numIS=100,
                                 inference.type="pre-multiply",
                                 improve.nomass.problem=TRUE, bits=1000,
                                 write.time=TRUE, verbose=TRUE)

    ## Append this to the problem
    return(list(resid.order=resid.order, intervals=intervals))
    ## Whenever there is a p-value fail, store whever is needed to rerun this.
    ## if(any(is.nan(results.rwbs))){

    ## } else {
    ##     return(pvs.rwbs)
    }
}


## Actually run the simulations.
nsim = 300 # 50
mc.cores = 5
results = mclapply(1:nsim,function(isim){
    printprogress(isim,nsim)
    return(onesim_rwbs(y.orig))
},mc.cores=mc.cores)
save(c("results", "resids"), file=file.path(outputdir, "artif-rwbs.Rdata"))


## Reading and seeing speed from file
read.time.from.file <- function(myfile){
    sort(readLines(myfile))
}

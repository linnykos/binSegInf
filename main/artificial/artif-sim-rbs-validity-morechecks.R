## Synopsis: We found that inference doesn't have correct type-I error. Why is
## this? Investigating what was going wrong..

## Synopsis: Now trying this with Gaussian residuals and assuming known noise
## level.
datadir = "../data"
outputdir = "../output"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))
library(microbenchmark)

## Verify from a smaller n=200 example, and with Gaussian noise
n = 100
lev = 3
sigma = 1
sigma.add=0.2
mn = fourjump(n=n, lev=lev)
sigma.add = sigma*0.2
resid = rnorm(n,0,sigma)

onesim_origgaussian_bootstrap_rbs <- function(y.orig, bits=1000, verbose=FALSE){

    ## Bootstrap data
    y = mn + bootstrap_sample(resid)

    ## Obtain inferences.
    how.close=5
    numIS=10
    min.num.things=30
    start.time = Sys.time()
    object = inference_bsFs(y=y, max.numSteps=15, consec=2, sigma=sigma,
                            postprocess= TRUE, locs=1:length(y), numIS= numIS,
                            min.num.things=min.num.things,
                            inference.type="pre-multiply", bits=bits,
                            sigma.add=sigma.add, verbose=verbose,
                            start.time=start.time, how.close=how.close,
                            retain.only.null.cases=TRUE,
                            mn=mn)
}

## Run actual simulations
bits = 2000
nsim = 5000
## mc.cores=3
mc.cores=3
start.time.overall = Sys.time()
results = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim,
                  lapsetime = round(difftime(Sys.time(), start.time.overall,
                                             units = "secs"), 2))
    ## Run a single result
   myresult = onesim_origgaussian_bootstrap_rbs(y.orig, bits=bits,  verbose=FALSE)
    ## if(any(myresult==1)) stop(msg=paste0("seed is ", isim))
}, mc.cores=mc.cores)

save(results, file=file.path(outputdir, "rbs-validity4.Rdata"))

# Synopsis: We found that inference doesn't have correct type-I error. Why is
# this? Investigating what was going wrong..
datadir = "../data"
outputdir = "../output"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))
library(microbenchmark)

onesim_gaussian_rbs <- function(y.orig, bits=1000, lev=1, verbose=FALSE, seed=NULL, old, max.numIS=2000, min.num.things=30){

    ## Verify from a smaller n=200 example, and with Gaussian noise
    n = 200
    lev = 3
    sigma = 1
    sigma.add=0.2
    mn = fourjump(n=n, lev=lev)
    if(!is.null(seed))set.seed(seed)
    y = mn + rnorm(n,0,sigma)

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
                            max.numIS=max.numIS,
                            mn=mn,
                            old=old)
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
   myresult = onesim_gaussian_rbs(y.orig, bits=bits,  verbose=FALSE)
    ## if(any(myresult==1)) stop(msg=paste0("seed is ", isim))
}, mc.cores=mc.cores)

save(results, file=file.path(outputdir, "rbs-validity4.Rdata"))
## validity 3 and 4 are identical except for bits=1000 goes to 2 000

## Also currently identifying when the ones occur.




## Load results and visualize things
## load(file=file.path(outputdir, "rbs-validity.Rdata"))
## load(file=file.path(outputdir, "rbs-valid

## Trying to fix isim=47
isim=47
myresult = onesim_gaussian_rbs(y.orig, bits=bits,  verbose=TRUE, seed=isim, old=FALSE, max.numIS=20000, min.num.things=1)
myresult = onesim_gaussian_rbs(y.orig, bits=4000,  verbose=TRUE, seed=isim, old=FALSE, max.numIS=20000, min.num.things=1)
myresult = onesim_gaussian_rbs(y.orig, bits=bits,  verbose=TRUE, seed=isim, old=TRUE)

## Some things take significantly longer.

## Would I get the same result with the old version of the code?



outputdir = "../output"
load(file=file.path(outputdir, "rbs-validity3.Rdata"))
pvs = (unlist(results[which(sapply(results,length)==1)]))
pvs = unlist(results[which(sapply(results, length)!=0)])
qqunif(pvs)

relevant.results = unlist(results[which(sapply(results, length)!=0)])
one.results = relevant.results[which(sapply(relevant.results, function(myresult) any(myresult==1)))]

one.results




## pvs.list = lapply(results, function(myresult){
##     myresult=myresult$results
##     pvs = lapply(myresult, function(myobj){myobj$pv})
## })

## pvs = Map(function(null,pvs){pvs[which(null==0)]}, null.list, pvs.list)
## qqunif(unlist(pvs))


## ## Load results after having done the null filtering.
## load(file=file.path(outputdir, "rbs-validity2.Rdata"))
## pvs.list = lapply(results, function(myresult){
##     myresult=myresult$results
##     pvs = lapply(myresult, function(myobj){myobj$pv})
## })
## qqunif(unlist(pvs))

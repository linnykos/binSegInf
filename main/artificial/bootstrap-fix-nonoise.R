## Synopsis: small-data examples to attempt to fix the problem of
## bootstrapping-noise procedure giving conservative p-value.
datadir = "../data"
outputdir = "../output"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))
library(microbenchmark)


n=40
subresid = resid.cleanmn[1:n]
## subresid = resid.cleanmn[1:n]
set.seed(0)
## subresid=rnorm(n,0,.1)
## sigma = sd(subresid)
## sigma=.1
sigmas = sapply(1:10000, function(isim){
        sd(bootstrap_sample(subresid))
})
sigma = mean(sigmas)

resid.cleanmn/sd(resid.clean.mn)


## Make a smaller problem
onesim_bootstrap_fix <- function(n=20, sigma=1,
                                sigma.add=sigma*0.2, bits=1000,
                                lev=1, verbose=FALSE, seed=NULL,
                                max.numIS=2000, min.num.things=30){
    ## Generate Gaussian noise
    mn = fourjump(n=n, lev=lev)
    y = mn + bootstrap_sample(subresid)

    ## Obtain inferences.
    postprocess = FALSE
    how.close = 1
    numIS = 10
    min.num.things = 30
    object = inference_bsFs(y=y, max.numSteps=15, consec=2, sigma=sigma,
                            postprocess= postprocess, locs=1:length(y), numIS= numIS,
                            min.num.things=min.num.things,
                            inference.type="pre-multiply", bits=bits,
                            sigma.add=sigma.add, verbose=verbose,
                            start.time=start.time, how.close=how.close,
                            max.numIS=max.numIS,
                            sim.options = list(retain.only.null.cases=TRUE,
                                               mn=mn,
                                               old=FALSE))
    return(object)
}



## Run the simulation
nsim = 500
mc.cores = 4
start.time = Sys.time()
results = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim, start.time=start.time)
    onesim_bootstrap_fix(n=20, sigma=sigma, sigma.add=0,
                        bits=5000, lev=0.1,
                        verbose=FALSE,
                        max.numIS=2000,
                        seed=isim,
                        min.num.things=30)}, mc.cores=mc.cores)
qqunif(unlist(results))

## Nonoise
## filename = "bootstrap-fix-small-gaus-nonoise.Rdata"
## filename = "bootstrap-fix-small-tdist-nonoise.Rdata"
## filename = "bootstrap-fix-small-tdist-nonoise-manydfs.Rdata"
## filename = "bootstrap-fix-small-orig-nonoise.Rdata"
## filename = "bootstrap-fix-small-orig-zeromean-nonoise.Rdata"
## filename = "bootstrap-fix-small-orig-zeromean-nonoise-cheatmore.Rdata"
## filename = "bootstrap-fix-small-orig-zeromean-nonoise-deflated.Rdata"
## filename = "bootstrap-fix-small-orig-nonoise-n200.Rdata"
## filename = "bootstrap-fix-small-orig-nonoise-n200-zeromean.Rdata"
## save(results, file=file.path(outputdir,filename))
## filename = "bootstrap-fix-small-tdist-nonoise-manydfs.Rdata"
## save(results.by.df, file=file.path(outputdir,filename))
## Parse the results
load(file=file.path(outputdir,filename))
makeplotfilename <- function(filename, ext= ".jpg"){ gsub(".Rdata", ext, filename )}
plotfilename = makeplotfilename(filename)
jpeg(file=file.path(outputdir, plotfilename), width=700,
     height=700)
## results1=results
## results.master = c(results1, results2)
## qqunif(unlist(results.master))
qqunif(unlist(results))
graphics.off()


## Data plot



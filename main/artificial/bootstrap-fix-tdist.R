## Synopsis: making a series of post-sel p-values from bootstrap, for n=20, but
## with varying degrees of freedom from df=n/4 to df=n.
datadir = "../data"
outputdir = "../output"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))
library(microbenchmark)

## Make a smaller problem
onesim_tdist_rbs <- function(n=200, sigma=1,
                                sigma.add=sigma*0.2, bits=1000,
                                lev=1, verbose=FALSE, seed=NULL,
                                max.numIS=2000, min.num.things=30,
                                dfreedom){
    ## Generate Gaussian noise
    mn = fourjump(n=n, lev=lev)
    if(!is.null(seed))set.seed(seed)
    y = mn + rt(n=n,df=dfreedom)
    sigma = sqrt(dfreedom/(dfreedom-2))

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

## Produce results
n = 20
dfs = c(5,10,20,40,80,160)
results.by.df <- lapply(dfs, function(dfreedom){

    printprogress(dfreedom, dfs, "degrees of freedom.",fill=TRUE)

    ## Run the simulation
    nsim = 2000
    mc.cores = 3
    start.time = Sys.time()
    results = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim, start.time=start.time)
        onesim_tdist_rbs(n=n, sigma=1, sigma.add=0,
                            bits=5000, lev=1,
                            verbose=FALSE,
                            max.numIS=2000,
                            min.num.things=30,
                            dfreedom=dfreedom)}, mc.cores=mc.cores)
})
filename = "bootstrap-tdist-nonoise-manydfs.Rdata"
save(results.by.df, file=file.path(outputdir, filename))

## Plot results
load(file=file.path(outputdir, filename))
makeplotfilename <- function(filename, ext= ".jpg"){ gsub(".Rdata", ext, filename )}
plotfilename = makeplotfilename(filename)
jpeg(file=file.path(outputdir, plotfilename), width=700,
     height=700)
par(mfrow=c(2,3))
for(ii in 1:6)qqunif(unlist(results.by.df[[ii]]))
graphics.off()


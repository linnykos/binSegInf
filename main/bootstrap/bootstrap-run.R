## Synopsis: script to run the bootstrap example. Run by ``Rscript
## ../bootstrap/bootstrap-run.R 1 2 3 4 5'', from the directory
## binSegInf/binSegInf

## Load the current package, helpers, and data.
load_all() 
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))
source(file=file.path("../main/bootstrap/bootstrap.R"))
outputdir = "~/Desktop/dummyoutput"


## Simulation settings
bits = 5000
mc.cores = 6
args = commandArgs(trailingOnly=TRUE)
facs = as.numeric(args)
nsims = seq(from=2000,to=5000, length=5)

## Over different sample sizes, collect results
for(fac in facs){
    printprogress(fac, facs, "factors", fill=TRUE)
    nsim = nsims[fac]
    start.time = Sys.time()
    results = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim,
                      lapsetime = round(difftime(Sys.time(), start.time,
                                                 units = "hours"), 2))
        myresult = onesim_rbs(y.orig, bits=bits, fac=fac, verbose=FALSE)
        return(myresult)
    }, mc.cores=mc.cores)
    
    ## Write to file
    facstring = paste0(unlist(strsplit(toString(fac), split='.', fixed=TRUE)), collapse="")
    filename = paste0("artif-rbs-fac-new-", facstring, ".Rdata")
    save(results, file=file.path(outputdir, filename))
}

## ## Check this before running it.
## microbenchmark::microbenchmark({
##     set.seed(1)
##     lev = 4
##     n = 100
##     mn = c(rep(0,n/2), rep(lev,n/2))
##     y = mn + rnorm(n,0,1)
##     sigma = 1; sigma.add = .2

##     start.time = Sys.time()
##     out = inference_bsFs(y=y, max.numSteps=15, consec=2,
##                          sigma=sigma, postprocess= TRUE,
##                          locs=1:length(y), numIS=10,
##                          min.num.things=30,
##                          inference.type="pre-multiply",
##                          bits=bits, sigma.add=sigma.add,
##                          verbose=verbose, start.time=start.time,
##                          how.close=5,
##                          mn=mn)
##     ## onesim_rbs(y.orig, bits=bits, fac=fac, verbose=FALSE)
## },times=5)

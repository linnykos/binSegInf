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
outputdir="../output"


## Simulation settings
bits = 5000
mc.cores = 6
args = commandArgs(trailingOnly=TRUE)
facs = as.numeric(args)
nsims = seq(from=2000,to=4000, length=5)

## Over different sample sizes, collect results
for(fac in facs){
    printprogress(fac, facs, "factors", fill=TRUE)
    nsim = nsims[fac]
    start.time = Sys.time()
    npart = 4
    for(part in 1:npart){
        results = mclapply(1:(nsim/npart), function(isim){
            printprogress(isim+(nsim/npart)*(part-1), (nsim),
                          lapsetime = round(difftime(Sys.time(), start.time,
                                                     units = "hours"), 2))
            myresult = onesim_rbs(y.orig, bits=bits, fac=fac, verbose=FALSE)
            return(myresult)
        }, mc.cores=mc.cores)
        
        ## Write to file
        ## facstring = paste0(unlist(strsplit(toString(fac), split='.', fixed=TRUE)), collapse="")
        filename = paste0("artif-rbs-fac-new-", fac, "-part-", part,".Rdata")
        cat("Saved", filename, fill=TRUE)
        save(results, file=file.path(outputdir, filename))
    }
}


## ## Timing this
## bits=5000
## reduced.time = microbenchmark::microbenchmark({
##     myresult = onesim_rbs(y.orig, bits=bits, fac=fac, verbose=TRUE, reduced=TRUE)
## }, times=5)

## nonreduced.time = microbenchmark::microbenchmark({
##     myresult = onesim_rbs(y.orig, bits=bits, fac=fac, verbose=TRUE, reduced=FALSE)
## }, times=5)



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


a=c(6, 
9, 
6, 
7, 
7, 
8, 
7, 
6, 
8, 
6, 
9, 
9, 
7, 
7, 
6, 
8, 
11, 
12, 
6, 
6, 
9, 
8, 
13, 
8, 
9, 
13, 
10, 
12, 
9, 
13, 
9, 
7, 
5, 
7, 
7, 
5, 
6, 
9, 
9, 
6, 
13, 
9, 
10, 
6, 
10, 
13, 
6, 
13, 
6, 
10, 
6, 
11, 
13, 
10, 
7, 
7, 
6, 
6, 
9, 
8, 
7, 
9, 
8, 
6, 
8, 
9, 
6, 
9, 
10, 
8, 
9, 
6, 
12, 
7, 
9, 
8, 
7, 
9, 
7, 
11)

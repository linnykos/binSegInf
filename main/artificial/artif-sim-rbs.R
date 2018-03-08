## Synopsis: The idea is to add bootstrapped residuals, then see the 1.
## conditional/unconditional powers of the tests done at the viscinities of the
## true guys, and 2. uniformity (or lack thereof) in the tests conducted
## regarding null locations.

# Data directory
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))

##' Multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1, n){
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}

## Simulation settings
onesim_rbs <- function(y.orig, bits=1000, fac=1, verbose=FALSE){

    ## Add bootstrapped residuals around a cleaned mean, with known sigma
    sigma = sd(y.orig[1:200]) * fac
    sigma.add = sigma*0.2
    y = newmn[-(1:200)] + bootstrap_sample(resid.cleanmn[-(1:200)]) * fac

    pvs.rbs = do_rbs_inference(y=y, max.numSteps=10, consec=2, sigma=sigma,
                               postprocess=TRUE, locs=1:length(y), numIS=100,
                               inference.type="pre-multiply", bits=bits, sigma.add=sigma.add,
                               write.time=TRUE, verbose=verbose)
    return(pvs.rbs)
}

## myresult = onesim_rbs(y.orig, bits=bits, fac=fac, verbose=TRUE)

## Factor
## fac = 2 ## 1, 1.5, 2, 2.5, 3
whichfac = 4
fac = c(1, 1.5, 2, 2.5, 3)[whichfac]
nsim = c(200,300,400,500,600)[whichfac]
## nsim = 200
bits = 3000
mc.cores = 7
results = list()
## for(isim in 1:nsim){
start.time = Sys.time()
results = mclapply(1:nsim, function(isim){

    printprogress(isim, nsim,
                  lapsetime = round(difftime(Sys.time(), start.time,
                                             units = "secs"), 2))
    ## Run a single result
    myresult = onesim_rbs(y.orig, bits=bits, fac=fac, verbose=FALSE)
}, mc.cores=mc.cores)

## Write to file
facstring = paste0(unlist(strsplit(toString(fac), split='.', fixed=TRUE)), collapse="")
filename = paste0("artif-rbs-fac-", facstring, "Rdata")
save(results, file.path(outputdir, filename))






## Reading and seeing speed from file
read.time.from.file <- function(myfile){
    sort(readLines(myfile))
}



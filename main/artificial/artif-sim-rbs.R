## Synopsis: The idea is to add bootstrapped residuals, then see the 1.
## conditional/unconditional powers of the tests done at the viscinities of the
## true guys, and 2. uniformity (or lack thereof) in the tests conducted
## regarding null locations.

# Data directory
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))
## Add bootstrapped residuals around a cleaned mean, with known sigma

## Simulation settings
onesim_rbs <- function(y.orig, bits=1000, fac=1, verbose=FALSE){

    ## Add bootstrapped residuals around a cleaned mean, with known sigma
    sigma = sd(y.orig[1:200]) * fac
    sigma.add = sigma*0.2
    y = newmn[-(1:200)] + bootstrap_sample(resid.cleanmn[-(1:200)]) * fac

    ## pvs.rbs = do_rbs_inference(y=y, max.numSteps=10, consec=2, sigma=sigma,
    ##                            postprocess=TRUE, locs=1:length(y), numIS=100,
    ##                            inference.type="pre-multiply", bits=bits, sigma.add=sigma.add,
    ##                            write.time=TRUE, verbose=verbose)
    ## return(pvs.rbs)

    how.close=5
    numIS=10
    min.num.things=30
    start.time = Sys.time()
    object = inference_bsFs(y=y, max.numSteps=15, consec=2,
                            sigma=sigma, postprocess= TRUE,
                            locs=1:length(y), numIS= numIS,
                            min.num.things=min.num.things,
                            inference.type="pre-multiply",
                            bits=bits, sigma.add=sigma.add,
                            verbose=verbose, start.time=start.time,
                            how.close=how.close
                            ## myloc=1859+c((-2):2)
                            ## myloc=923+c((-2):2)
                            )
}

## myresult = onesim_rbs(y.orig, bits=bits, fac=fac, verbose=TRUE)

## Factor
fac=1
myresult = onesim_rbs(y.orig, bits=bits, fac=fac, verbose=TRUE)


nsim = c(200,250,300,350,400)[fac]
## nsim = 200
bits = 1000
mc.cores = 6
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
filename = paste0("artif-rbs-fac-new-", facstring, ".Rdata")
save(results, file=file.path(outputdir, filename))



## Attempt to see intermediate results; to be applied in next iteration.
facstring = paste0(unlist(strsplit(toString(fac), split='.', fixed=TRUE)), collapse="")
filename = paste0("artif-rbs-fac-", facstring, "-chunk.Rdata")
nchunk = 5
nsim.chunk = nsim/nchunk
results.list = list()
for(ichunk in 1:nchunk){
    chunkresult = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim,
                      lapsetime = round(difftime(Sys.time(), start.time,
                                                 units = "secs"), 2))
        ## Run a single result
        myresult = onesim_rbs(y.orig, bits=bits, fac=fac, verbose=FALSE)
    }, mc.cores=mc.cores)
    results.list[[ichunk]] = chunkresult
    save(results.list, file=file.path(outputdir, filename))
}

## Combining results
save(results.list, file=file.path(outputdir, filename))

## Synopsis: Bootstrap + IC stopping does not work. why?
datadir = "../data"
outputdir = "../output"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))

## ## Residuals
subresid = resid.cleanmn[1:n]
## subresid = rnorm(n,0,1)
sigma_orig= sd(subresid)
nsim=10000
a = sapply(1:nsim, function(isim){
    y = rep(0,n) + bootstrap_sample(subresid)
    sd(y)
})
sigma=mean(a)


## Make a smaller problem
onesim <- function(n=200, sigma=1,
                   sigma.add=sigma*0.2, bits=1000, verbose=FALSE, seed=NULL,
                   max.numIS=2000, min.num.things=30){
    mn = rep(0,n)
    subresid = resid.cleanmn[1:n]
    ## sigma = sd(subresid)
    ## mn = fourjump(n=n,lev=.1)
    set.seed(seed)
    bootstrapped.resid = bootstrap_sample(subresid)
    ## sigma = sd(bootstrapped.resid)
    y = mn + bootstrapped.resid

    ## Obtain inferences.
    postprocess = FALSE
    how.close = 1
    numIS = 10
    min.num.things = 30
    object = inference_bsFs(y=y, max.numSteps=13, consec=2, sigma=sigma,
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
nsim = 300
mc.cores = 4
start.time = Sys.time()
## sigma=mean(a)
results = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim, start.time=start.time)
    onesim(n=20, sigma=sigma, sigma.add=0,
                        bits=5000,
                        verbose=FALSE,
                        max.numIS=2000, seed=isim+nsim,
                        min.num.things=30)}, mc.cores=mc.cores)

qqunif(unlist(results))

## Synopsis: small-data examples to attempt to fix the problem of
## bootstrapping-noise procedure giving conservative p-value.
datadir = "../data"
outputdir = "../output"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))
library(microbenchmark)

## Make a smaller problem
onesim_bootstrap_rbs <- function(len=200, bits=1000,  verbose=FALSE, seed=NULL,
                                 max.numIS=2000, min.num.things=30){
    if(!is.null(seed))set.seed(seed)
    inds = sample(length(resid.cleanmn), len)
    mn = rep(0, len)
    y = mn + bootstrap_sample(resid.cleanmn[inds])
    sigma.add=0

    ## Obtain inferences.
    postprocess = FALSE
    how.close = 1
    numIS = 10
    min.num.things = 30
    verbose=FALSE
    object = inference_bsFs(y=y, max.numSteps=1, consec=2, sigma=sigma,icstop=FALSE,
                            postprocess= postprocess, locs=1:length(y), numIS= numIS,
                            min.num.things=min.num.things,
                            inference.type="pre-multiply", bits=bits,
                            sigma.add=sigma.add, verbose=verbose,
                            start.time=start.time, how.close=how.close,
                            max.numIS=max.numIS,
                            sim.options = list(retain.only.null.cases=FALSE,
                                               mn=len,
                                               old=FALSE))
    return(object)
}

mc.cores = 4
nsim = 1000
sigma = sd(resid.cleanmn)
start.time = Sys.time()
results = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim, start.time=start.time)
    onesim_bootstrap_rbs(len=200, bits=5000,  verbose=FALSE,
                        max.numIS=2000, min.num.things=30)
},mc.cores=4)
filename = "bootstrap-rbs-full.jpg"
save(results, file=file.path(outputdir,filename))

sigma
qqunif(unlist(results))




## Run the simulation
## results.list = lapply(c(1, 1.1, 1.2, 1.3),function(fac){
results.list = lapply(c(1, 2, 3, 4, 5)+5, function(seed){

## set.seed(0)
## n=200
## n=400
n=100
## set.seed(6)
set.seed(seed)
subresid = resid.cleanmn[sample(length(resid.cleanmn), n)]
## subresid = resid.cleanmn[sample(length(resid.cleanmn), n)]
    fac=1
sigma = sd(subresid)*fac
plot(subresid)
nsim = 500
mc.cores = 4
start.time = Sys.time()
results = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim, start.time=start.time)
    onesim_bootstrap_rbs(n=n, sigma=sigma, sigma.add=0,#0.2*sigma,
                        bits=5000, lev=.1, subresid=subresid,
                        verbose=FALSE,
                        max.numIS=2000,
                        seed=isim,
                        min.num.things=30,
                        fac=fac)}
  , mc.cores=mc.cores)
})
## results
qqunif(unlist(results))

filename = "bootstrap-fix-n200-orig-nonoise-diff.jpg"
filename = "bootstrap-fix-n200-orig-nonoise-differentseeds.jpg"
filename = "bootstrap-orig.jpg"
jpeg(file=file.path(outputdir, filename), width=700,height=700)
## qqunif(unlist(results))
## graphics.off()
par(mfrow=c(2,3))
for(ifac in 1:5){
    ## fac = c(1, 1.1, 1.2, 1.3)[ifac]
    result = results.list[[ifac]]
    qqunif(unlist(result))
    set.seed(ifac)
    subresid = resid.cleanmn[sample(length(resid.cleanmn), n)]
    ## plot(subresid, ylim=range(resid.cleanmn))
    ## abline(h=mean(subresid))
    title(main=round(sd(subresid),3))
}
graphics.off()


## filename = "bootstrap-orig-nonoise-n200.Rdata" ## Hydra3
filename = "bootstrap-orig-yesnoise.Rdata" ## Hydra3
## filename = "bootstrap-orig-nonoise-zeromean.Rdata" ## Hydra5
## filename = "bootstrap-orig-yesnoise.Rdata" ## Hydra4
save(results, file=file.path(outputdir,filename))


filename = "bootstrap-orig-nonoise-n200.Rdata" ## Hydra3
filename = "bootstrap-orig-yesnoise-n200.Rdata" ## Hydra3
outputdir="../output"
load(file=file.path(outputdir,filename))
makeplotfilename <- function(filename, ext= ".jpg"){ gsub(".Rdata", ext, filename )}
plotfilename = makeplotfilename(filename)
jpeg(file=file.path(outputdir, plotfilename), width=700,
     height=700)
qqunif(unlist(results))
graphics.off()


subresid = resid.cleanmn[1:n]
subresid = subresid-mean(subresid)
sigma_orig= sd(subresid)
nsim=10000
sigmas = sapply(1:nsim, function(isim){
    sd(bootstrap_sample(subresid))
})
mean(sigmas)
plot(density(sigmas))
abline(v=sigma_orig, lwd=3)
abline(v=mean(sigmas), lwd=3, col='blue')



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

## Additive noise
## filename = "bootstrap-fix-small-gaus-yesnoise.Rdata"
## filename = "bootstrap-fix-small-tdist-yesnoise.Rdata"
## filename = "bootstrap-fix-n200-orig-yesnoise.Rdata"
## filename = "bootstrap-fix-small-orig-zeromean-yesnoise-deflated.Rdata"
## filename = "bootstrap-fix-small-orig-zeromean-yesnoise.Rdata"
## save(results, file=file.path(outputdir,filename))

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

## Seeing if the standard deviation of the bootstrapped residuals actually
## converges to that of the original sample.
n=200
subresid = resid.cleanmn[1:n]
subresid = rnorm(n,0,1)##subresid - mean(subresid) ## centering
sd.orig = sd(subresid)
sd.orig2 = mysd(subresid)
nsim=10000
sds = sapply(1:nsim, function(isim){
    bs = bootstrap_sample(subresid)
    sigma = sd(bs)
    return(sigma)
})
plot(density(sds))
abline(v=sd.orig,lwd=3)
## sigma.add = sigma*0.2

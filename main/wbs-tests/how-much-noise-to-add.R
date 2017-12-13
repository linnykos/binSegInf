##' Synopsis: explore the effect of additive noise on the recovery and power of
##' binseg.
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
## nsim=12
## n=20
dosim_additive_noise <- function(sigma.add=sigma.add, n=20, nsim, lev=1){
    visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
    results = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim)
        pvs = dosim_compare(type="sbs.rand", n=n, lev=lev, numIntervals=n,
                            sigma.add=sigma.add, numIS=100, meanfun=fourjump,
                            visc=visc.fourjump, numSteps=4, bits=1000)
    }, mc.cores=4)
    return(results)
}

ngrain = 5
sigma.add.list = seq(from=0, to=3, length=ngrain)
nsims = seq(from=1000, to=100, len=ngrain)/10
results.list = list()
## for(isigma in 1:length(sigma.add.list)){
results = lapply(1:ngrain, function(igrain){
    printprogress(igrain, ngrain)
    sigma.add = sigma.add.list[igrain]
    dosim_additive_noise(sigma.add=sigma.add, n=20, nsim=nsims[[igrain]])
})


## Visualize conditional powers
par(mfrow=c(2,6))
for(ii in 1:10){
    results = results.list[[ii]]
    qqunif(unlist(sapply(results, function(myresult)myresult[,"pvs"])))
}

## Visualize the detection rate
for(ii in 1:10){
    results = results.list[[ii]]
    ## qqunif(unlist(sapply(results, function(myresult)myresult[,"pvs"])))
    sum(sapply(results, nrow))
}

##' Synopsis: explore the effect of additive noise on the recovery and power of
##' binseg.
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
## nsim=12
## n=20
dosim_additive_noise <- function(sigma.add=sigma.add, n=20, nsim, lev=1, mc.cores=4){

    ## Sample settings
    sigma.add=1
    nsim=10
    n = 20
    lev=1
    nsim=12

    visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
    pvs.list = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim)
        isim=1
        set.seed(isim)
        pvs = dosim_compare(type="sbs.rand", n=n, lev=lev, numIntervals=n,
                            sigma.add=sigma.add, numIS=100, meanfun=fourjump,
                            visc=visc.fourjump, numSteps=4, bits=1000,
                            max.numIS=10000)
    }, mc.cores=mc.cores)
    return(pvs.list)
}

ngrain = 5
sigma.add.list = seq(from=0, to=2, length=ngrain)
nsims = rep(2000, ngrain)

## sigma.add.list = seq(from=0, to=3, length=ngrain)
## nsims = seq(from=1000, to=100, len=ngrain)/10
## ## results.list = list()
## ## for(isigma in 1:length(sigma.add.list)){

results = lapply(1:ngrain, function(igrain){
    printprogress(igrain, ngrain, type="sigma grain")
    cat(fill=TRUE)
    sigma.add = sigma.add.list[igrain]
    dosim_additive_noise(sigma.add=sigma.add, n=200, nsim=nsims[[igrain]],
                         mc.cores=8)
})

save(sigma.add.list, nsims, results,
     file=file.path(outputdir, "how-much-noise-to-add.Rdata"))

load(file=file.path(outputdir, "how-much-noise-to-add.Rdata"))

## Format and plot pvalues
pvs.list = lapply(results, function(myresults){
    pvs = sapply(myresults, function(myresult)myresult[,"pvs"])
    pvs = unlist(pvs)
    pvs = pvs[!is.na(pvs)]
    pvs = pvs[pvs!=1]
})
qqunif(pvs.list, cols=1:5)
legend("bottomright", col=1:5, lty=rep(2,5), lwd=rep(2,5),legend=sigma.add.list  )


## Visualize the detection rate
detections = c()
sums=c()
for(ii in 1:10){
    myresults = results[[ii]]
    sums[ii] = sum(sapply(myresults, nrow))
    ## detections[ii] = (length(unlist(locs.list)) / (4*length(locs.list)))
}
## plot(detections, type='l')
plot(sums, type='l')

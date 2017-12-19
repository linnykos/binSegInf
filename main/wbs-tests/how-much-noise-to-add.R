##' Synopsis: explore the effect of additive noise on the recovery and power of
##' binseg.
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
## nsim=12
## n=20
dosim_additive_noise <- function(sigma.add=sigma.add, n=20, nsim, lev=1, mc.cores=4){
    visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
    pvs.list = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim)
        pvs = dosim_compare(type="sbs.rand", n=n, lev=lev, numIntervals=n,
                            sigma.add=sigma.add, numIS=100, meanfun=fourjump,
                            visc=visc.fourjump, numSteps=4, bits=1000)
    }, mc.cores=mc.cores)
    return(pvs.list)
}

ngrain = 5
sigma.add.list = seq(from=0, to=2, length=ngrain)
nsims = rep(2000, ngrain)
results = lapply(1:ngrain, function(igrain){
    printprogress(igrain, ngrain, type="sigma grain")
    cat(fill=TRUE)
    sigma.add = sigma.add.list[igrain]
    dosim_additive_noise(sigma.add=sigma.add, n=200, nsim=nsims[[igrain]],
                         mc.cores=8)
})


## Visualize conditional powers
par(mfrow=c(2,6))
for(ii in 1:10){
    myresult = results[[ii]]
    qqunif(unlist(sapply(myresult, function(mything)mything[,"pvs"])))
}


## Visualize the detection rate
detections = c()
for(ii in 1:5){
    myresult = results[[ii]]
    ## qqunif(unlist(sapply(myresult, function(mything)mything[,"pvs"])
    locs.list = sapply(myresult, function(mything)mything[,"locs"])
    detections[ii] = (length(unlist(locs.list)) / (4*length(locs.list)))
}
plot(detections, type='l')

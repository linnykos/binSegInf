##' Synopsis: explore the effect of additive noise on the recovery and power of
##' binseg.
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
## nsim=12
## n=20
dosim_additive_noise <- function(sigma.add=sigma.add, n=20, nsim, lev=1, mc.cores=4){

    ## ## Timing a single thing
    ## sigma.add=1
    ## nsim=10
    ## n = 20
    ## lev=1
    ## nsim=12
    ## microbenchmark::microbenchmark({
    ## pvs = dosim_compare(type="sbs.rand", n=n, lev=lev, numIntervals=n,
    ##                     sigma.add=sigma.add, numIS=100, meanfun=fourjump,
    ##                     visc=visc.fourjump, numSteps=4, bits=1000,
    ##                     max.numIS=10000)
    ## print(pvs)
    ## }, times = 5)

    ## visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
    visc.fourjump = 1:n
    pvs.list = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim)
        ## set.seed(isim)
        pvs = dosim_compare(type="sbs.rand", n=n, lev=lev, numIntervals=n,
                            sigma.add=sigma.add, numIS=100, meanfun=fourjump,
                            visc=visc.fourjump, numSteps=4, bits=3000,
                            max.numIS=20000) ## This takes either really small or up to 10 minutes.
    }, mc.cores=mc.cores)
    return(pvs.list)
}

ngrain = 10
sigma.add.list = seq(from=0, to=2, length=ngrain)
nsim = 2000
nsims = rep(nsim, ngrain)

## results = list()
## for(igrain in 1:ngrain){
## for(igrain in 5:7){
for(igrain in 8:10){
    printprogress(igrain, ngrain, type="sigma grain")
    cat(fill=TRUE)
    sigma.add = sigma.add.list[igrain]
    results[[igrain]] = dosim_additive_noise(sigma.add=sigma.add, n=200, nsim=nsims[[igrain]],
                         mc.cores=7)
    save(sigma.add.list, nsims, results,
         file=file.path(outputdir, "how-much-noise-to-add-moresim8910.Rdata"))
    ## save(sigma.add.list, nsims, results,
         ## file=file.path(outputdir, "how-much-noise-to-add-moresim567.Rdata"))
}

## load(file=file.path(outputdir, "how-much-noise-to-add.Rdata"))##ngrain=5,nsim=300
load(file=file.path(outputdir, "how-much-noise-to-add-moresim.Rdata"))

## Format and plot conditional pvalues
pvs.list = lapply(results, function(myresults){
    pvs = sapply(myresults, function(myresult)myresult[,"pvs"])
    pvs = unlist(pvs)
    pvs = pvs[!is.na(pvs)]
    pvs = pvs[pvs!=1]
})
qqunif(pvs.list, cols=1:5)
legend("bottomright", col=1:5, lty=rep(2,5), lwd=rep(2,5),legend=sigma.add.list  )

cond.powers = c()
detections = c()
for(igrain in 1:ngrain){

    myresults = results[[igrain]]

    ## Get conditional power
    myverdicts = sapply(myresults, function(myresult){
        numtest = nrow(myresult)
        return(myresult[,"pvs"] < 0.05/numtest )
    })
    sumtests = length(unlist(myverdicts))
    cond.powers[igrain] = sum(unlist(myverdicts), na.rm=TRUE)/sumtests

    ## Get detection
    detections[igrain] = sum(sapply(myresults, nrow))/(nsims[igrain]*4)
}

plot(detections~sigma.add.list, lwd=2, ylim=c(0, 1), type='l')
lines(cond.powers~sigma.add.list, lwd=2, col='red')
legend("bottomleft", col=c("black", "red"), lwd=c(2,2), legend = c("detection", "conditional power"))

## To be fair, I need to make unconditional powers as well.




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

pvs.list


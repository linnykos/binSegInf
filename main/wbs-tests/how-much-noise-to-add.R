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

## Calculate the unconditional powers
cond.pows = c()
for(ii in 1:5){
    myresults = results[[ii]]
    numtests = sapply(myresults, nrow)
    rejs = sapply(myresults, function(mymat){
        pvs = mymat[,"pvs"]
        pvs[is.nan(pvs)] <- 0
        rejections = sum(pvs < 0.05/4)
    })
    uncond.pows[ii] = sum(rejs)/sum(numtests)
}

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




## Visualize unconditional power
n = 200
visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-2,-1,0,1,2)))
cond.pows = detections = c()
for(ii in 1:5){
    myresults = results[[ii]]
    totalnumtests = sum(sapply(myresults, nrow))
    rejs = sapply(myresults, function(mymat){
        pvs = mymat[,"pvs"]
        locs = mymat[,"locs"]
        pvs = pvs[which(locs %in% visc.fourjump)]
        pvs[is.nan(pvs)] <- 0
        rejections = sum(pvs < 0.05/4)
    })
    numtests = sapply(myresults, function(mymat){
        locs = mymat[,"locs"]
        sum(locs %in% visc.fourjump)
    })
    print(summary((numtests)))
    cond.pows[ii] = sum(rejs)/sum(numtests)
    detections[ii] = sum(numtests)/totalnumtests
}

## Visualize the cond/uncond powers and detection, all together
col.cond = "black"
col.uncond = "red"
col.detection = "black"
lty.detection = 2
w=5; h=5
pdf(file=file.path(outputdir, "powers-how-much-noise-to-add.pdf"), width=w, height=h)
mar = c(4.5,4.5,0.5,0.5)
par(mar=mar)
plot(NA, ylim=c(0,1), type='l', lwd=2, xlim = range(sigma.add.list), axes=FALSE,
     ylab = "Unconditional power", xlab = bquote(Additive~noise~sigma[add]))
axis(1); axis(2)
lines(uncond.pows~sigma.add.list[1:5], type='o', col=col.uncond, lwd=2)
lines(cond.pows~sigma.add.list[1:5], type='o', col=col.cond, lwd=2)
lines(detections~sigma.add.list[1:5], type = 'o', col=col.detection, lty=lty.detection, lwd=2)
graphics.off()


## Make the QQ plots
mycols = RColorBrewer::brewer.pal(5, "Set2")
w=5; h=5
mar = c(4.5,4.5,0.5,0.5)
pdf(file=file.path(outputdir, "qqplot-how-much-noise-to-add.pdf"), width=w, height=h)
par(mar=mar)
qqunif(pvs.list, cols=mycols[1:5])
legend("bottomright", col=mycols[1:5], lty=rep(2,5), lwd=rep(2,5),legend=round(sigma.add.list[1:5],3))
graphics.off()




## ## Visualize the detection rate
## detections = c()
## sums=c()
## ## for(ii in 1:10){
## for(ii in 1:5){
##     myresults = results[[ii]]
##     sums[ii] = sum(sapply(myresults, nrow))
##     ## detections[ii] = (length(unlist(locs.list)) / (4*length(locs.list)))
##     ## detections[ii] = (length(unlist(locs.list)) / (4*length(locs.list)))
##     myresults[[1]]
## }
## ## plot(detections, type='l')
## plot(sums, type='l')

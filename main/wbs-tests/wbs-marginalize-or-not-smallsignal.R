## Synopsis: Want to see if marginalization of WBS inference is worth it, in
## terms of (conditional) power! Rerunning experiments for smaller signal

outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
dosim_compare_rwbs_to_wbs <- function(n=20, numIntervals=n, nsim, lev=1,
                                      mc.cores=4, type = c("wbs.nonrand", "wbs.rand"),
                                      visc=1:n){

    ## visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
    type = match.arg(type)
    start.time = Sys.time()
    all.times = rep(NA,nsim)
    pvs.list = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim,
                      lapsetime = round(difftime(Sys.time(), start.time,
                                            units = "secs"), 2))
        pvs = dosim_compare(type=type,
                            n=n, lev=lev, numIntervals=numIntervals,
                            numIS=100, meanfun=fourjump,
                            visc=visc, numSteps=4,
                            max.numIS=20000)
        return(pvs)
    }, mc.cores=mc.cores)
    return(pvs.list)
}





## Running this code on a server
jj = 2
n = 200
numIntervals = n
mc.cores = 8
whichlev.list = list(1:2, 3:4, 5:6)
whichlev = whichlev.list[[jj]]
levs = c(0.1,0.2,0.3,0.4,0.5, 0.75)[whichlev]
nsims = seq(from=3000,to=1500,length=6)[whichlev]
mc.cores = 8
visc.firstjump = n/5 + c(-1,0,1)
visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
results.wbs = results.rwbs = list()
filename = paste0("wbs-vs-rwbs-lev-", paste0(whichlev, collapse=""), ".Rdata")
print(filename)
for(ilev in 1:length(levs)){
    lev = levs[ilev]
    nsim = nsims[ilev]

    ## Non-Marginalized WBS
    results.wbs[[ilev]] = dosim_compare_rwbs_to_wbs(numIntervals=numIntervals, n=n, lev=lev,
                                           nsim=nsim, mc.cores=mc.cores, type = "wbs.nonrand",
                                           visc=visc.firstjump)
    save(results.rwbs, results.wbs, file=file.path(outputdir, filename))

    ## Non-Marginalized WBS
    results.rwbs[[ilev]] = dosim_compare_rwbs_to_wbs(numIntervals=numIntervals, n=n, lev=lev,
                                            nsim=nsim, mc.cores=mc.cores, type = "wbs.rand",
                                            visc=visc.firstjump)
    save(results.rwbs, results.wbs, file=file.path(outputdir, filename))
}


## Load and aggregate the results
whichlev.list = list(1:2, 3:4, 5:6)
results.rwbs.master = results.wbs.master = list()
for(jj in 1:3){
    whichlev = whichlev.list[[jj]]
    filename = paste0("wbs-vs-rwbs-lev-", paste0(whichlev, collapse=""), ".Rdata")
    outputdir="../output"
    load(file=file.path(outputdir, filename))

    results.wbs.master[whichlev] = results.wbs
    if(jj==3){
        results.rwbs.master[[5]] = results.rwbs[[1]]
    } else {
        results.rwbs.master[whichlev] = results.rwbs
    }
}

## results.wbs.master[whichlev] = results.wbs
## results.wbs.master[[5]] = results.wbs
## results.rwbs.master[[5]] = results.rwbs[[1]]
levs = c(0.1,0.2,0.3,0.4,0.5, 0.75)[1:5]

results.rwbs = results.rwbs.master
results.wbs = results.wbs.master

## Reformating to extract /just/ the p-values.
pvs.list.wbs = pvs.list.rwbs = list()
for(ilev in 1:length(levs)){
    myresult = results.wbs[[ilev]]
    pvs = unlist(sapply(myresult, function(a)a$pvs))
    pvs.list.wbs[[ilev]] = pvs[which(!is.na(pvs))]

    myresult = results.rwbs[[ilev]]
    pvs = unlist(sapply(myresult, function(a)a$pvs))
    pvs.list.rwbs[[ilev]] = pvs[which(!is.na(pvs))]
}


## Pairwise qqplot
pdf(file=file.path("~/Desktop", "smallsignal.pdf"))
cols = RColorBrewer::brewer.pal(length(pvs.list.wbs), "Set1")
par(mfrow=c(2,3))
for(ilev in 1:length(levs)){
    qqunif(list(wbs=pvs.list.wbs[[ilev]], rwbs=pvs.list.rwbs[[ilev]]), cols=cols)
    title(main=paste0("lev=", levs[ilev]))
}
graphics.off()
n=200

lev=0.5
mn = fourjump(lev,n)
sigma=1
y = mn + rnorm(n, 0, sigma)
plot(y, pch=16, col='grey50', main = expression(delta==0.5))
lines(mn, lwd=3, col='red')
abline(v=n/5, col='lightgrey', lty=2)
graphics.off()



## Just rerunning 0.03 (many cores)
npart = 4
n = 200
numIntervals = n
lev=0.3
mc.cores = 8
visc.firstjump = n/5 + c(-1,0,1)
nsim = 1500
ipart = 4
filename = paste0("wbs-vs-rwbs-repeat-part-", ipart, ".Rdata")
print(filename)
results.wbs.parts = results.rwbs.parts = list()

## Non-Marginalized WBS
a = dosim_compare_rwbs_to_wbs(numIntervals=numIntervals, n=n, lev=lev,
                              nsim=nsim, mc.cores=mc.cores, type = "wbs.nonrand",
                              visc=visc.firstjump)
results.wbs.parts[[ipart]] = a

## Non-Marginalized WBS
nsim=500 ## 1500
b = dosim_compare_rwbs_to_wbs(numIntervals=numIntervals, n=n, lev=lev,
                                                       nsim=nsim, mc.cores=mc.cores, type = "wbs.rand",
                                                 visc=visc.firstjump)
results.rwbs.parts[[ipart]] = b
save(results.rwbs.parts, results.wbs.parts, file=file.path(outputdir, filename))



## Load and aggregate results
results.wbs.parts.master = results.rwbs.parts.master = list()
for(ipart in 1:npart){
    ipart=1
    filename = paste0("wbs-vs-rwbs-repeat-part-", ipart, ".Rdata")
    load(file=file.path(outputdir, filename))
    results.wbs.parts.master[[ipart]] = results.wbs.parts[[ipart]]
    results.rwbs.parts.master[[ipart]] = results.rwbs.parts[[ipart]]
}

## Combine previous results ase well.
filename = paste0("wbs-vs-rwbs-lev-34.Rdata")
outputdir="../output"
load(file=file.path(outputdir, filename))
results.wbs.parts.master[[5]] = results.wbs[[1]]
results.rwbs.parts.master[[5]] = results.rwbs[[1]]


pv.wbs = unlist(sapply(results.wbs.parts.master, function(mything){
    sapply(mything, function(a){a$pvs})
}))

pv.rwbs = unlist(sapply(results.rwbs.parts.master, function(mything){
    sapply(mything, function(a){a$pvs})
}))
pv.wbs = pv.wbs[which(!is.na(pv.wbs))]
pv.rwbs = pv.rwbs[which(!is.na(pv.rwbs))]

w=h=5
pdf(file=file.path("~/Desktop", "zeropointthree.pdf"), width=w, height=h)
qqunif(list(wbs=pv.wbs, rwbs=pv.rwbs), cols=cols)
graphics.off()
length(pv.wbs)
length(pv.rwbs)

## Power difference (doesn't show that much because not a lot of rejections at
## all in such a small signal regime.)
sum(pv.wbs<0.05/4)
sum(pv.rwbs<0.05/4)


## ## Separate calculation of detection abilities in small signal, for $n=200$
## source("../main/wbs-tests/sim-helpers.R")
## sigma = 1
## n = 200
## meanfun = fourjump
## nsim = 2000
## numIntervals = n
## numSteps = 4
## visc.firstjump = n/5 + c(-1,0,1)
## levs = c(0.1,0.2,0.3,0.4,0.5, 0.75)[whichlev]
## rwbs.detections = sapply(levs, function(lev){
##     printprogress(lev, levs, "levels")
##     cat(fill=TRUE)
##     retains = mclapply(1:nsim, function(isim){
##         printprogress(isim, nsim)
##         mn = meanfun(lev=lev,n=n)
##         y = mn + rnorm(n, 0, sigma)
##         cumsum.y = cumsum(y)
##         inference.type = "pre-multiply"
##         improve.nomass.problem = TRUE
##         g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
##         retain = which(abs(g$cp) %in% visc.firstjump)
##         return(length(retain))
##     }, mc.cores=8)
##     return(sum(unlist(retains))/nsim)
## })
## names(rwbs.detections) = levs

## ## Save as file
## filename="wbs-detection-smallsignal.Rdata"
## save(rwbs.detections, file=file.path(outputdir, filename))

## ## Load this





## Synopsis: See how stable a given p-value is (by fixing data and intervals, and only
## repeating the importance sampling)

outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")

## Generate data once
n = 200
lev = 1 ## lev=2
## lev=2
sigma = 1
mn = fourjump(lev=lev,n=n)
y = mn + rnorm(n, 0, sigma)
plot(y)
points(mn,lwd=3)
cumsum.y = cumsum(y)
inference.type = "pre-multiply"
improve.nomass.problem = TRUE
numIS = 100
max.numIS=2000
numIS=100
bits=3000
min.num.things = 30

## Apply WBS
numSteps = 4
set.seed(1)
intervals = intervals(n=n, numIntervals=n)
g = wildBinSeg_fixedSteps(y, intervals=intervals, numSteps=numSteps)
poly.wbs = polyhedra(obj=g$gamma, u=g$u)
vlist <- make_all_segment_contrasts(g)
locs = as.numeric(names(vlist))

## Get the p-values
mc.cores = 8
v = vlist[[1]]
cumsum.v = cumsum(v)
cumsum.y = cumsum(y)
nsim = 100
min.num.things.list = c(10, 30,50,70,90,100)
infolists = list()
for(ithing in 1:length(min.num.things.list)){
    printprogress(ithing, length(min.num.things.list), "min num things")
    cat(fill=TRUE)
    min.num.things = min.num.things.list[ithing]
    start.time = Sys.time()
    infolist = list()
    for(isim in 1:nsim){
        printprogress(isim,nsim, "total simulations",
                      lapsetime = round(difftime(Sys.time(), start.time,
                                                 units = "secs"), 2))
        infolist[[isim]] = randomize_wbsfs(v=v, winning.wbs.obj=g, sigma=sigma,
                                           numIS=numIS,
                                           inference.type=inference.type,
                                           cumsum.y=cumsum.y,cumsum.v=cumsum.v,
                                           improve.nomass.problem=improve.nomass.problem,
                                           bits=bits, max.numIS=max.numIS,
                                           return.more.things=TRUE, min.num.things=min.num.things,
                                           mc.cores=mc.cores, verbose=FALSE)
    }
    infolists[[ithing]] = infolist
    ## Save results
    ## filename = "how-stable-is-rwbs.Rdata"
    filename = "how-stable-is-rwbs-harder-problem.Rdata"
    save(list=ls(), file=file.path(outputdir, filename))
}




## Load results
filename = "how-stable-is-rwbs.Rdata"
filename = "how-stable-is-rwbs-harder-problem.Rdata"
load(file=file.path(outputdir, filename))

## Plot the six histograms of the p-values with different stop criteria for
## importance sampling.
w=10;h=5;
## pdf(file=file.path(outputdir, "how-stable-is-rwbs-hist.pdf"), width=w, height=h)
pdf(file=file.path(outputdir, "how-stable-is-rwbs-hist-harder.pdf"), width=w, height=h)
par(mfrow=c(2,3))
par(mar=c(2,1,1,1))
for(ii in 1:length(infolists)){
    infolist = infolists[[ii]]
    pvs = sapply(infolist, function(a)a$pv)
    ## if(ii==1)pvs=pvs[pvs<1E-10] ## for the non-harder, original case
    hist(pvs, xlim=c(0,1), main = min.num.things.list[ii], breaks=50) ## xlim = c(0,1)##xlim=c(0,4E-12)
    ## plot(density(pvs),xlim=c(0,1))
    print(range(pvs))
}
graphics.off()

## Plot the spreads (fitted standard deviations).
spreads = c()
for(ii in 1:length(infolists)){
    infolist = infolists[[ii]]
    spreads[ii] = (sd(sapply(infolist, function(a)a$pv)))
}
w=h=5;
## pdf(file=file.path(outputdir, "how-stable-is-rwbs-sd.pdf"), width=w, height=h)
pdf(file=file.path(outputdir, "how-stable-is-rwbs-sd-harder.pdf"), width=w, height=h)
plot((spreads)~min.num.things.list[1:length(infolists)], type = 'o', lwd=2,
     pch=16, cex=2, xlim= range(min.num.things.list),
     main = "Sd of p-values, \n by min. # of contributing IS draws")
graphics.off()
## plot(spreads~min.num.things.list[1:length(infolists)], type = 'o', lwd=2, pch=16, cex=2,
##      main = "Sd of p-values, by min. # of nonzero things")


## Plot data and contrast, for progress-slides/binSegInf/slides.tex
w=h=5;
pdf(file=file.path(outputdir, "how-stable-is-rwbs-data-plot-harder.pdf"), width=w, height=h)
## pdf(file=file.path(outputdir, "how-stable-is-rwbs-data-plot.pdf"), width=w, height=h)
par(mar=c(2,1,1,1))
plot(y, pch=16, col = 'grey50', axes=FALSE)
axis(2);axis(1)
lines(mn, col='red', lwd=3)
v2=v
v2[v==0]=NA
points(v2*10, col='blue', cex=2, pch=16)
graphics.off()


## How many replicates are required, for each case?
w=h=5;
## pdf(file=file.path(outputdir, "how-stable-is-rwbs-numdraws.pdf"), width=w, height=h)
pdf(file=file.path(outputdir, "how-stable-is-rwbs-numdraws-harder.pdf"), width=w, height=h)
mean.draws = sd.draws = c()
for(ii in 1:length(infolists)){
    infolist = infolists[[ii]]
    mean.draws[ii] = mean(sapply(infolist, function(a)a$numIS.cumulative))
    sd.draws[ii] = sd(sapply(infolist, function(a)a$numIS.cumulative))
}
mymat = cbind(mean.draws+sd.draws,mean.draws, mean.draws-sd.draws)
matplot(NA, ylim = range(as.numeric(mymat)), col='white', xlim = range(min.num.things.list), axes=FALSE, xlab = "Min # of valid IS draws", ylab = "total # of draws that were run")
axis(1);axis(2)
for(irow in 1:nrow(mymat)){
    myrow = mymat[irow,]
    print(myrow)
    lowtick = myrow[1]
    midtick =  myrow[2]
    uptick = myrow[3]
    x.to.plot = min.num.things.list[irow]
    segments(x0=x.to.plot,x1=x.to.plot,y0=lowtick,y1=uptick)
    points(x=x.to.plot,y=lowtick, pch="-")
    points(x=x.to.plot,y=midtick, pch=16)
    points(x=x.to.plot,y=uptick, pch="-")
}
title(main = "How many IS draws are needed for \n certain # of valid IS draws")
graphics.off()

## The reason this spread is small is because it is a strong signal anyway??
## I'll be trying on a more challenging problem (lev=1)
## Also, there is a lot of more information that we /could/ use, like how long it took!! i.e. IS efficiency. Plot this

## I should probably do this on a larger n, like n=2000 or n=200

## Also, I should do it on more challenging problems, which are probably
## problems in which the polyhedra are more restrictive.


## Settings to explroe
## (1) one-jump, one-step case
## (2) one-jump, two-step case
## (3) four-jump, four-step case


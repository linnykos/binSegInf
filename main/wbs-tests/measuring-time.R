## Synopsis: measure time for randomized WBS simulations.
library(microbenchmark )

## Computations by sample size
ns = c(50, 100, 200, 300,  500)
slow.times.by.n = list()
for(ii in 1:length(ns)){
    print(ii)
    n = ns[ii]
    numIntervals = n

    ## Generate some data
    lev=0
    mn = onejump(lev,n)
    set.seed(0)
    y = mn + rnorm(n, 0, sigma)
    cumsum.y = cumsum(y)
    numsteps=1

    g.orig = wildBinSeg_fixedSteps(y, numIntervals=numIntervals,
                                   numSteps= numSteps,
                                   inference.type="pre-multiply")
    vlist <- make_all_segment_contrasts(g.orig)
    numIS=100
    v = vlist[[1]]
    cumsum.v = cumsum(v)
    slow.times.by.n[[ii]] = microbenchmark({
        randomize_wbsfs(v, winning.wbs.obj=g.orig, numIS = 100, sigma=sigma,
                        comprehensive=FALSE, inference.type="pre-multiply",
                        cumsum.y=cumsum.y, cumsum.v=cumsum.v)
    }, times=3)
}


myplot <- function(times.by.n, ns, ylim = c(0,20000), add=FALSE, col='black',...){
    getmedtime <- function(mytime){(mytime$time[3])/1000000000}
    medtimes= sapply(times.by.n, getmedtime)
    ns.squared = ns^2
    dat = data.frame(medtimes, ns, ns.squared)
    g = lm(medtimes ~ ns + ns.squared, data=dat)
    newdat = data.frame(ns = 1:2000, ns.squared = (1:2000)^2)
    proj.times = predict(g, newdata=newdat)
    if(add) myplotter = points
    if(!add) myplotter = plot
    myplotter(medtimes ~ ns, xlim=c(0,2000), ylim = ylim,col=col,...)
    new.ns=  1:2000
    lines(proj.times~new.ns,col=col)
}

pdf("~/Desktop/times-after-functionmagic.pdf", width=5, height=5)
myplot(times.by.n, ns, pch=16, col='black',cex=1.5)
myplot(slow.times.by.n, ns, add=TRUE, pch=16, col='red',cex=1.5)
legend("topleft", legend = c("before change", "after change"), pch=c(16,16), col=c("black", "red"),cex=c(1.5,1.5))
graphics.off()


pdf("~/Desktop/times-after-functionmagic2.pdf", width=5, height=5)
myplot(times.by.n, ns, pch=16, col='black',cex=1.5, ylim=c(0,600))
graphics.off()



## Is the p-value distribution any more/less stable in high $n$? I think it will
## be more unstable, so we need to grow p-value; I'd be (pleasantly surprised)
## if it weren't.



ns = c(5,10,20)
distr.by.n = list()
for(ii in 1:length(ns)){
    print(ii)
    n = ns[ii]
    numIntervals = n

    ## Generate some data
    lev=0
    mn = onejump(lev,n)
    set.seed(0)
    y = mn + rnorm(n, 0, sigma)
    cumsum.y = cumsum(y)
    numsteps=1
    g.orig = wildBinSeg_fixedSteps(y, numIntervals=numIntervals,
                                   numSteps= numSteps,
                                   inference.type="pre-multiply")
    vlist <- make_all_segment_contrasts(g.orig)
    numIS=50
    v = vlist[[1]]
    cumsum.v = cumsum(v)
    nsim=100
    distr.by.n[[ii]] = mclapply(1:nsim, function(isim){
        printprogress(isim,nsim)
        randomize_wbsfs(v, winning.wbs.obj=g.orig, numIS = 100, sigma=sigma,
                        comprehensive=FALSE, inference.type="pre-multiply",
                        cumsum.y=cumsum.y, cumsum.v=cumsum.v)
    }, mc.cores=3)
}

## Load
n=60
lev=0
sigma = 1
## Generate some data
mn = onejump(lev,n)
set.seed(0)
y = mn + rnorm(n, 0, sigma)
numsteps=1

## Fit WBS
numIntervals = n
set.seed(1)
g.orig = wildBinSeg_fixedSteps(y, numIntervals=numIntervals,
                               numSteps= numSteps,
                               inference.type="rows")
poly = polyhedra(obj=g.orig$gamma, u=g.orig$u)
vlist <- make_all_segment_contrasts(g.orig)
numIS=100
v = vlist[[1]]
poly.pval2(y=y, poly=poly, v=v, sigma=sigma)$pv

cumsum.y = cumsum(y)
cumsum.v = cumsum(v)

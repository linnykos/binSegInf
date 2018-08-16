##Synopsis: Among replicates, the stoptimes were sometimes very large, but 90%
##of the time, they were at most 7.
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
sigma = sd(y.orig[1:200])
nsim = 100
stoptimes = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim)
    set.seed(isim)
    max.numSteps = 20
    consec = 2
    y = y.orig[201:length(y.orig)]
    numIntervals=length(y.orig)
    g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=max.numSteps,
                              inference.type='none')
    ic_obj = get_ic(g$cp, y, consec=consec, sigma=sigma, type="bic")
    return(ic_obj$stoptime)
}, mc.cores=3)
table(stoptimes)


## Is it an egregious problem that stoptimes can get really big? <- If we can do
## various things like decluttering to remedy it, then we might not worry about
## this. Otherwise, it would be great to be able to stabilize this. By taking
## the median/average stoptimes over subsampled (or equidistant harvested)
## partitions of the data.

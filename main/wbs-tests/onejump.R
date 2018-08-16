## Synopsis: one-jump inference examples
source("../main/justin/plot-helpers.R")

## Helper to get p-values
onejumpsim = function(lev, n, nsim){
    cat("lev=", lev, fill=TRUE)
    numSteps = 3
    sigma = 1

    results = mclapply(1:nsim,function(isim){
        printprogress(isim, nsim)

        ## Generate some data
        mn = c(rep(0,n/2), rep(lev,n/2))
        y = mn + rnorm(n, 0, sigma)

        ## Fit WBS
        numIntervals = n
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
        poly = polyhedra(obj=g$gamma, u=g$u)
        vlist <- make_all_segment_contrasts(g)
        pvs = sapply(vlist, function(v){
            return(poly.pval2(y=y, poly=poly, v=v, sigma=sigma)$pv)
        })
        names(pvs) = g$cp*g$cp.sign

        ## Also store other information
        null.true = sapply(vlist, function(v){
            return(v%*%mn == 0)
        })

        return(list(pvs=pvs, null.true=null.true))
    },mc.cores=4)
    cat(fill=TRUE)

    pvs = unlist(lapply(results, function(a)a[["pvs"]]))
    truths = unlist(lapply(results, function(a)a[["null.true"]]))
    return(list(pvs=pvs, truths=truths))
}

## Actually run simulations
levs = c(0,1,2,3)
results = lapply(levs, onejumpsim, n=60, nsim=10)
pvs = results$pvs
truths = results$truths

## Save results
outputdir = "../output"
filename = "onejump.Rdata"
save(list=c("levs","pvs","locs","truths"), file=filename)

## Load results
load(file=file.path(outputdir,filename))

## Plot + Export results
pdf("nonrand-wbs.pdf")

## Process results
visc = (n/2+((-3):3))
visc = (n/2+((-1):1))
visc = (n/2)
locs = lapply(pvs, function(mypvs) as.numeric(names(mypvs)))
visc.pvs <- Map(function(mypv, myloc){
    mypv[myloc %in% (n/2+((-3):3))]},pvs, locs)
names(cond.pvs) = paste("jump-size=",levs)

qqunif(cond.pvs,cols=1:4)
maketitle(n=n, other="Three steps")

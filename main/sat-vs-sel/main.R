source('../main/sat-vs-sel/sim-helper.R')

## Generate a toy dataset
n = 12
twojump <- function(lev,n){c(rep(0,n/3),rep(lev,n/3), rep(0,n/3))}
lev=1
sigma=1
mn = twojump(lev,n)

nsim = 100
results = mclapply(1:nsim, function(isim){
    printprogress(isim,nsim)
    y = mn + rnorm(n, 0, sigma)

    ## Fit WBS
    numSteps=2
    ## g = binSeg_fixedSteps(y, numSteps=numSteps)
    g = wildBinSeg_fixedSteps(y, numSteps=numSteps, numIntervals=n)
    poly = polyhedra(obj=g$gamma, u=g$u)
    vlist <- make_all_segment_contrasts(g)

    ## Saturated model p-values
    sat.pvs = sapply(vlist, function(v){
        poly.pval2(y=y,v=v,poly=poly,sigma=sigma)$pv
    })

    ## Selected model p-values
    cps = abs(as.numeric(names(vlist)))
    cp.rest.list = sapply(cps, function(mycp)(cps)[cps!=mycp])
    covariance = diag(sigma,n)
    sel.pvs = unlist(Map(function(v,cp.rest){
        sel.pvals(y=y, v=v, poly=poly, sigma=sigma, cp.rest = cp.rest, covariance=covariance,ngen=100000)
    }, vlist, cp.rest.list))

    return(list(sat.pvs=sat.pvs, sel.pvs=sel.pvs))
}, mc.cores=4)

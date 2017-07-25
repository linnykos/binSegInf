## Synopsis: We want to see uniform p-values when
## Simulation settings.
numSteps = 1
nsim = 100
nsim.is = 100
lev = 0
## pvec = rep(NA,nsim)

## Generate data
cat(fill=TRUE)
## for(isim in 1:nsim){
pvec = mclapply(1:nsim, function(isim){
    cat('\r', isim, "out of", nsim)
    my.mn <- meanfun(lev,n)
    y = (my.mn + stats::rnorm(n,0,sigma))

    ## Do fused lasso inference
    D = makeDmat(n,type='tf',ord=0)
    obj <- genlassoinf::dualpathSvd2(y,D=D,maxsteps=numSteps,approx=TRUE)
    contrasts <- make_all_segment_contrasts(obj)
    pvec = pvec.plain = setNames(rep(NA,length(obj$cp)), obj$cp)
    poly <- polyhedra(obj$Gobj.naive$G, obj$Gobj.naive$u)

    ## pv <- poly.pval2(y=y, poly=poly, v=contrasts[[ii]], sigma=sigma,
    ##                              reduce=reduce)$pv
    pv = randomized_genlasso_pv(y=y, v=contrasts[[ii]],
                                       sigma=sigma, numSteps=numSteps,
                                       numIntervals=numIntervals,
                                       nsim.is=nsim.is, bits=100,
                                       reduce=reduce, augment=augment)
    return(pv)

}, mc.cores = 3)

hist(unlist(pvec),xlim=c(0,1))

## Do fused lasso inference.

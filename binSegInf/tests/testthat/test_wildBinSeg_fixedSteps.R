test_that("WBS.FT gives uniform p-values",{

    nsim=3000
    pvs = lapply(1:nsim, function(isim){

        ## Generate some data
        n = 10
        lev = 0
        mn = c(rep(0,n/2), rep(lev,n/2))
        sigma = 1
        y = mn + rnorm(n, 0, sigma)

        ## Fit WBS
        numIntervals = 3
        numSteps = 3
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
        poly = polyhedra(obj=g$gamma, u=g$u)
        vlist <- make_all_segment_contrasts(g)
        return(sapply(vlist, function(v){
            return(poly.pval2(y=y, poly=poly, v=v, sigma=sigma)$pv)
        }))
    })
    expect_true(ks.test(unlist(pvs), punif)$p.value > 0.05)
})

test_that("WBS has power when signal is present.",{
    lev = 5
    nsim = 1000
    n = 10
    numSteps = 3
    sigma = 1

    pvs = mclapply(1:nsim,function(isim){

        ## Generate some data
        mn = c(rep(0,n/2), rep(lev,n/2))
        set.seed(isim)
        y = mn + rnorm(n, 0, sigma)

        ## Fit WBS
        numIntervals = 30
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
        poly = polyhedra(obj=g$gamma, u=g$u)
        vlist <- make_all_segment_contrasts(g)
        pvs = sapply(vlist, function(v){
            return(poly.pval2(y=y, poly=poly, v=v, sigma=sigma)$pv)
        })
        names(pvs) = g$cp*g$cp.sign
        return(pvs)
    },mc.cores=4)
    pvs = unlist(pvs)
    pvs = unique(pvs) ## Sometimes ties occur numerically. Not sure why.
    expect_true(ks.test(unlist(pvs), punif)$p.value < 0.01)
})


test_that("WBS polyhedron is exact",{

    ## Simulation settings
    n = 10
    lev = 0
    mn = c(rep(0,n/2), rep(lev,n/2))
    sigma = 1
    numIntervals = 30

    nsim=1000
    a = mclapply(1:nsim,function(isim){

        ## Generate original data.
        y.orig = mn + rnorm(n, 0, sigma)
        obj.orig = intervals(numIntervals=numIntervals, n=n)

        ## Fit WBS
        numSteps = 3
        g.orig = wildBinSeg_fixedSteps(y.orig, numIntervals=numIntervals, numSteps=numSteps, intervals=obj.orig)
        poly.orig = polyhedra(obj=g.orig$gamma, u=g.orig$u)

        ## Add noise and see if exact correspondence
        sigma.add = .2
        y.new = y.orig + rnorm(n,0,sigma.add)
        g.new = wildBinSeg_fixedSteps(y.new, numIntervals=numIntervals, numSteps=numSteps, intervals=obj.orig)

        ## Helper function to see if results match up, at every step.
        wbs_results_match <- function(obj.orig, obj.new){
            return(all.equal(obj.orig$results[,c("max.s", "max.b", "max.e", "max.sign")],
                             obj.new$results[,c("max.s", "max.b", "max.e", "max.sign")]))
        }
        ## See if the results are exactly the same.
        if(wbs_results_match(g.orig, g.new)==TRUE){
            expect_true(contained(poly.orig, y.new))
        }
    }, mc.cores=4)
})

test_that("Randomized wbs inference has uniform null p-values",{

    ## Simulation settings
    n = 4
    lev = 0
    mn = c(rep(0,n/2), rep(lev,n/2))
    sigma = 1
    numIntervals = 3
    numSteps=1

    ## Get randomized p-values
    nsim = 1000
    pvs = mclapply(1:nsim, function(isim){
        printprogress(isim,nsim)
        set.seed(isim)
        y = mn + rnorm(n, 0, sigma)
        pv = suppressWarnings(randomize_wbsft(y=y, numSteps=numSteps, numIntervals=numIntervals, numIS = 20,
                                            comprehensive=TRUE)[[1]])
    }, mc.cores=4)
    ## qqunif(unlist(pvs))
})

test_that("WBS output matched wbs::wbs()",{
print('not written yet')
})

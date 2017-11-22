context("Test IC wrapper functions")

test_that("Object returned from ic_to_polyhedra() is a valid polyhedra", {

    ## Example
    set.seed(1)
    y = c(runif(10),runif(10)+5,runif(10))
    sigma = .5
    cp = c(10,20,15,5)
    consec=1
    obj <-  get_ic(cp=cp, y=y, sigma=sigma,consec=consec, maxsteps=3)
    poly <- ic_to_poly(obj)
    expect_true(is_valid.polyhedra(poly))
})






## test_that("IC wrapper works properly", {

##     ## Example
##     set.seed(1)
##     y = c(rnorm(10),rnorm(10)+5,rnorm(10))
##     sigma = 1
##     obj = binSeg_fixedThresh(y,1)
##     cp = c(10,20,15,5)
##     consec=1
##     stp <- ic_wrapper(obj,y,sigma=sigma)
##     obj <-  get_ic(cp=cp, y=y, sigma=sigma,consec=consec, maxsteps=3)
##     poly <- ic_to_poly(obj)
##     expect_true(is_valid.polyhedra(poly))
## })



test_that("IC minimization works properly", {

    nsim=20
    lev = 0
    n = 10
    meanfun=onejump
    sigma=1
    numIntervals=n
    mn = meanfun(lev,n)
    set.seed(21)
    y = mn + rnorm(n, 0, sigma)
    cumsum.y = cumsum(y)

    ## Fit initial WBS for a generous number of steps
    numSteps=10
    g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps,
                              inference.type='rows')

    ## Get ic object of two steps
    consec=2
    ic_obj = get_ic(g$cp, g$y, consec=consec, sigma=sigma, type="bic")
    ic_poly = ic_obj$poly
    ic_poly$gamma%*%y > ic_poly$u
    expect_true(all(ic_poly$gamma %*% g$y >= ic_poly$u))

    ## Generate new y's
    nsim=1000
    sigmadd=1
    for(isim in 1:nsim){
        ynew = mn + rnorm(n, 0, sigma)
        ic_obj_new = get_ic(g$cp, ynew, consec=consec, sigma=sigma, type="bic")
        if(all(ic_poly$gamma %*% ynew >= ic_poly$u)){
            expect_equal(ic_obj_new$stoptime, ic_obj$stoptime)
            expect_equal(ic_obj_new$seqdirs, ic_obj$seqdirs)
        } ## else {
        ##     if(!is.na(ic_obj_new$stoptime) & ic_obj_new$stoptime==2){
        ##         print(ic_obj_new$stoptime)
        ##         print(ic_obj_new$seqdir)
        ##     }
        ## }
    }
})

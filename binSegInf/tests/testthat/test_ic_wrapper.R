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






test_that("IC wrapper works properly", {

    ## Example
    set.seed(1)
    y = c(rnorm(10),rnorm(10)+5,rnorm(10))
    sigma = 1 
    obj = binSeg_fixedThresh(y,1)

    cp = c(10,20,15,5) 

    consec=1
    obj <-  get_ic(cp=cp, y=y, sigma=sigma,consec=consec, maxsteps=3)
    poly <- ic_to_poly(obj)
    expect_true(is_valid.polyhedra(poly))
})

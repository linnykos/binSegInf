context("Test the functions binseg.OOO()")
require(wbs)


test_that("No NAs are returned from the cusum() function.", {
    n = 10
    sigma=.1
    y = rnorm(n,0,sigma)
    cusums = cusum(s=1,e=n,b=5,y=y,contrast = TRUE)
    expect_equal(any(is.na(cusums)), FALSE)
})

test_that("Changepoints match with those from wbs::sbs()", {

    ## Generate data
    n = 30 
    sigma=1
    lev1=0
    lev2=3
    set.seed(0)
    y = c(rep(lev1,n/2),rep(lev2,n/2)) + rnorm(n,0,sigma)
    thresh =  1.5
    
    ## Our changepoint set using SBS with fixed threshold
    blist = binseg.by.thresh(s = 1, e = n, j = 0, k = 1, thresh = thresh, y = y)$blist
    our.cpt=sort(blist[,3])
    
    ## Their changepoint set
    a = sbs(y)
    a.cpt = changepoints(a,th=thresh)
    their.cpt = sort(a.cpt$cpt.th[[1]])

    ## Match them.
    expect_equal(our.cpt, their.cpt)
})

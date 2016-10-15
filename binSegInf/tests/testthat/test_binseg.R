context("Test the basic standard binary segmentation function binseg()")
require(wbs)


## test_that("Changepoints match with those from wbs::sbs()", {

##     ## Generate data
##     n = 30 
##     sigma=1
##     lev1=0
##     lev2=3
##     set.seed(0)
##     y = c(rep(lev1,n/2),rep(lev2,n/2)) + rnorm(n,0,sigma)
##     thresh =  1.5
    
##     ## Run SBS with fixed threshold
##     slist = elist = blist = Blist = zlist = Zlist = matrix(NA, nrow = n, ncol = 2^8)
##     binseg(s = 1, e = n, j = 0, k = 1, thresh = thresh, y = y, n = n)
    
##     ## Our changepoint set
##     our.cpt=sort(as.numeric(na.omit(as.numeric(trim(blist)))))
    
##     ## Their changepoint set
##     a = sbs(y)
##     a.cpt = changepoints(a,th=thresh)
##     their.cpt = sort(a.cpt$cpt.th[[1]])

##     ## Match them.
##     expect_equal(our.cpt, their.cpt)
## })

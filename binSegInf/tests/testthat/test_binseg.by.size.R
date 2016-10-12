## context("Test binseg.by.size()")


## test_that("Works with /full/ number of steps", {
##     set.seed(0)
##     n=1000
##     for(n in c(10,100,1000,5000)){
##         set.seed(n)
##         y = rnorm(n,0,1)
##         a = binseg(y,n)
##         bs.output = list(blist = matrix(1,2,3,4))
##         expect_error(make.v(test.b, bs.output))
##     }
## })

## test_that("Works on null signal", {
##     set.seed(0)
##     y = rnorm(10,0,1)
##     binseg(y,4)
##     bs.output = list(blist = matrix(1,2,3,4))
##    expect_error(make.v(test.b, bs.output))
## })


## test_that("Works on null signal", {
##     set.seed(0)
##     y = rnorm(10,0,1)
##     binseg(y,4)
##     bs.output = list(blist = matrix(1,2,3,4))
##    expect_error(make.v(test.b, bs.output))
## })



## test_that("Works on relly jagged/noisy y", {
##     set.seed(1)
##     y = rnorm(10,0,10)
##     bs.output = list(blist = matrix(1,2,3,4))
##    expect_error(make.v(test.b, bs.output))
## })


## test_that("Doesn't work if there are NAs in y", {
##     test.b = 5 
##     bs.output = list(blist = matrix(1,2,3,4))
##    expect_error(make.v(test.b, bs.output))
## })



## test_that("Resulting contrast doesn't make a negative v^Ty (for one-sided testing)",{
##    test.b = 5
##    bs.output = list(blist = matrix(5),
##                     zlist = matrix(1),
##                     y = c(rep(5,5), rep(-5,5)))
##    expect_error(make.v(test.b, bs.output))
## })

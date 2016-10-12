context("Test contrast-forming function make.v()")

test_that("Coordinate that is tested is contained in blist", {
    test.b = 5 
    bs.output = list(blist = matrix(1,2,3,4))
   expect_error(make.v(test.b, bs.output))
})


test_that("Resulting contrast doesn't make a negative v^Ty (for one-sided testing)",{
   test.b = 5
   bs.output = list(blist = matrix(5),
                    zlist = matrix(1),
                    y = c(rep(5,5), rep(-5,5)))
   expect_error(make.v(test.b, bs.output))
})

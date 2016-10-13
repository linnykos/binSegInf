context("Compare binseg() and binseg.by.size() output.")

n = 12
set.seed(0)
y = c(rnorm(n/3,0,.5),rnorm(n/3,3,.5), rnorm(n/3,5,.5))

## Run it two ways
a = binseg.by.size(y, n-1, verbose=FALSE)
b = binseg.by.thresh(y, 0,verbose=FALSE)

## Compare the output
test_that("The changepoint sets are the same", {
    expect_equal(sort(a$B),
                 sort(as.numeric(trim(b$blist))))
})

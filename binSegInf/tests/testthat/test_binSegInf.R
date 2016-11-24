context("Compare binseg() and binseg.by.size() output.")

n = 12
set.seed(0)
y = c(rnorm(n/3,0,.5),rnorm(n/3,3,.5), rnorm(n/3,5,.5))

## Run /full/ standard binseg two ways
b = binseg.by.thresh(y, 0,verbose=FALSE)
a = binseg.by.size(y, n-1, verbose=FALSE)

## Compare the output
test_that("The changepoint sets are the same from two methods", {
    expect_equal(sort(a$B),
                 sort(as.numeric(trim(b$cp))))
})


## Compare the inference output.
test_that("(A sort of weak test that) Inference under a null signal is valid, by K-S test against uniformity.", {
    nsim = 100
    n = 12
    theta = rep(0,n)
    sigma = 1
    mysimulate = function(n,sigma){
        y = theta + rnorm(n, 0, sigma)
        b = binseg.by.size(y,3,verbose=FALSE)
        v = make.v(b$B[1],b$B,b$Z,n)
        return(poly.pval(y,b$G,b$u,v,sigma)$pv)
    }
    ps = replicate(nsim,mysimulate(n,sigma))
    expect_equal((ks.test(ps, "punif")$p.value > 0.05), TRUE) 
})

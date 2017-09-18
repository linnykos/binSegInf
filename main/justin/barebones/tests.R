test_that("Maximizer is found correctly.", {
    y = c(rep(0,5),rep(5,5))
    expect_equal(estim(y),5)
})

test_that("Polyhedron is exactly collect.", {
    n=10
    y = rep(0,n) + rnorm(n,0,0.1)
    cp = estim(y)
    nsim=1000
    for(isim in 1:nsim){
        ynew = y + rnorm(n,0,0.1)
        mypoly = polyhedron(cp, n)
        if(all(mypoly %*% ynew >= 0)){
            expect_equal(estim(ynew), cp)
        } else {
            expect_false(isTRUE(estim(ynew)== cp))
        }
    }
})

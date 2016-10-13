context("Test where.jk()")

test_that("It works on a manually constructed matrix A", {
    A = matrix(NA, nrow=40, ncol=3)
    set.seed(0)
    A[,1] = rep(1:10,each=4)
    A[,2] = as.numeric(replicate(10, sample(1:10,4,replace=F)))
    A[,3] = sample(1:100,40)
    expect_equal(where.jk(5,8,A),19)
})

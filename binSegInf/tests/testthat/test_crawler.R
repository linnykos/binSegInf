context("Test crawler().")
test_that("The crawler finds the right end indices for which myvec==mygoal", {
    myvec = c(1:4,rep(5,5),6:10)
    n = 12
    myloc=7
    mygoal=5
    expect_equal(crawler(myvec,myloc,mygoal), c(5:9))
})

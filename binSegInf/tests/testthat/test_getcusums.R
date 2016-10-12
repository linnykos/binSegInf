context("Test getcusum()")


test_that("has consistently signed output (in terms of $sign, $cusum, $abs.cusum", {
    y = rnorm(n,0,1)
    obj = getcusums(1,length(y),y)
    expect_equal(max(abs(obj$allcusums)), obj$cusum)
    expect_equal((abs(obj$allcusums)),obj$signs*obj$allcusums)
})

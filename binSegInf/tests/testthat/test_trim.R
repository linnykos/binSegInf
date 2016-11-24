context("Test trimming function trim()")

test_that("Trims all-NA matrices gracefully (by returning NULL)", {
    mymat = matrix(NA,ncol=10,nrow=10)
    expect_equal(trim(mymat), NULL)
})

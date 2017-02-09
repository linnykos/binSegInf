context("Test wbs.R functions")

## test .computeJumpIdx

test_that("Intervals are collected correctly", {
    
    ## Test settings
    s = 1
    e = 100
    numIntervals = 100

    ## See if it always gives the right number of elements.
    replicate(100, { expect_equal(length(generate_intervals(s, e,
        numIntervals)), numIntervals) })

    ## See if it produces NULL elements.
    replicate(100, {
        intervals = generate_intervals(s, e, numIntervals)
        expect_false(any(unlist(lapply(intervals, is.null))))
    })

})

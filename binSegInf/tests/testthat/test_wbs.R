context("Test wbs.R and wbs_helper.R functions.")

test_that("Intervals are collected correctly", {
    
    ## Test settings
    s = 1
    e = 100
    numIntervals = 1000

    ## See if it always gives the right number of elements.
    expect_equal(length(generate_intervals(s, e, numIntervals)$se),
                 numIntervals)

    ## See if it produces NULL elements.
    intervals = generate_intervals(s, e, numIntervals)
    expect_false(any(unlist(lapply(intervals, is.null))))
})


test_that("get_which_qualify() is correctly functioning", {

    set.seed(0)
    s=1
    e=60
    intervals = generate_intervals(s,e,10)
    m = which(get_which_qualify(15,45,intervals))
    expect_equal(sort(m), c(1,2,9))
})


test_that("make_se() is correctly functioning",{

    ## Test setting.
    s=1
    e=60
    y = rnorm(rep(0,4,each=30),0,1)
    intervals = generate_intervals(s,e,10)
    m = which(get_which_qualify(s,e,intervals))
    semat = .make_se_mat(m,intervals,y)

    ## Separately obtain all cusums.
    all.s = lapply(intervals$se[m], function(se) se[1])
    all.e = lapply(intervals$se[m], function(se) se[2])
    all.max.cusums.manual = unlist(Map(cusum, all.s, semat[,"b"], all.e, rep(list(y),length(m))))

    ## See
    expect_equal(semat[,"maxcusum"], all.max.cusums.manual)

})


test_that(".make_se_mat() is correctly functioning", {

    ## Test setting
    s=1
    e=60
    y = rnorm(rep(0,4,each=30),0,1)
    intervals = generate_intervals(s,e,10)
    m = which(get_which_qualify(s,e,intervals))
    semat = .make_se_mat(m,intervals,y)

    ## Separately obtain all cusums
    all.s = lapply(intervals$se[m], function(se) se[1])
    all.e = lapply(intervals$se[m], function(se) se[2])
    all.max.cusums.manual = unlist(Map(cusum, all.s, semat[,"b"], all.e, rep(list(y),length(m))))

    ## See if the maximizing is correctly done, with respect to internal and
    ## externally created max cusums.
    expect_equal(as.numeric(abs(semat[which(semat[,"maxhere"]==1), "maxcusum"])),
                 max(abs(semat[,"maxcusum"])))
    expect_equal(abs(all.max.cusums.manual[which(semat[,"maxhere"]==1)]),
                 max(abs(semat[,"maxcusum"])))
})

test_that("wbs() doesn't produce environment |env| whose tree element env$tree has null elements",{


    ## Test setting
    s=1
    e=60
    y = rnorm(rep(0,4,each=30),0,1)
    env = wbs(y,1,return.env=TRUE)

    ## See if the tree has any empty (NULL) elements
    expect_true(all(lapply(env$tree, length)>0))

})


test_that("wbs() gives the correct changepoint in a strong-signal case.", {

    ## Test setting
    thresh=10
    seed=0
    set.seed(seed)
    y = rep(c(0,10),each=30) + rnorm(60,0,1)
    numInterval=100

    ## Run WBS
    output = wbs(y, thresh, numInterval,seed=seed)

    ## See if single upward changepoint is detected
    expect_equal(as.numeric(output$cp), 30)
    expect_equal(as.numeric(output$cp.sign), +1)
    
}

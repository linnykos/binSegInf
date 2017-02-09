context("Testing the function rid_jk_list_nicely().")

test_that("Under normal settings, it deletes elements and makes elements scoot back, so that all nonzero elements of Tcurr come first and are consecutive, and all NULLs are trailing and consecutive."{
    Tcurr = list(c(1,1),c(1,2),NULL,NULL)
    Ecurr = Scurr = cplist(10)
    Ecurr = add(Ecurr, 1,1,7)
    Scurr = add(Scurr, 1,1,8)
    Ecurr = add(Ecurr, 1,2,11)
    Scurr = add(Scurr, 1,2,13)
    expect_equal(rid.jk.nicely(Tcurr,Ecurr,Scurr), list(c(1,2),NULL,NULL,NULL))
})


test_that("When it gets rid of the single element from a Tcurr object with one element, then it returns a list full of NULLs, of the same length.", {

    Tcurr = list(c(3,3),NULL,NULL)
    Ecurr = cplist(10)
    Ecurr = add(Ecurr, 3,3,7)
    Scurr = cplist(10)
    Scurr = add(Scurr, 3,3,8)
    expect_equal(rid.jk.nicely(Tcurr,Ecurr,Scurr), list(NULL,NULL,NULL))
})


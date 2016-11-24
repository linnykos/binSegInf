context("Test that cplist methods work correctly")

## Make a cplist
methods(class="cplist")
a=cplist(4)
a = add(a,1,1,1/3)
a = add(a,2,3,2/3)
a = add(a,3,5,5/3)


test_that("You can find valid existing indices the cplist.", {
    expect_equal(where_jk(a,1,1),1)
    expect_equal(where_jk(a,2,3),2)
    expect_equal(where_jk(a,3,5),3)
})

test_that("You can't find something that isn't in the cplist.", {
    expect_equal(where_jk(a,4,4),NA)
})



test_that("You aren't able to add (j,k,val) if (j,k) already exists in the matrix.", {
    expect_error(add(a,1,1,10/3))
})


test_that("extract.cplist() works properly",{
    expect_equal(extract(a,1,1), 1/3)
})


test_that("extract.cplist() doesn't find something that's not in there",{
   expect_error(extract(a,1,3))
   expect_error(extract(a,10,1))
})

test_that("trim.cplist() works properly",{
   trim.cplist(a)
   trim(a)
   expect_error(extract(a,10,1))
})

test_that("mydeparse() separates a string of the form 10:32 properly",{
    expect_error(mydeparse("10::32"))
    expect_equal(mydeparse("10:32"), 10:32)
})

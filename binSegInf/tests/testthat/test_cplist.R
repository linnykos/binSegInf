context("Test that cplist methods work correctly")

## Make a cplist
methods(class="cplist")
a=cplist(4)
a = add.cplist(a,1,1,1/3)
a = add.cplist(a,2,3,2/3)
a = add.cplist(a,3,5,5/3)


test_that("You can find valid existing indices the cplist.", {
    expect_equal(where_jk.cplist(a,1,1),1)
    expect_equal(where_jk.cplist(a,2,3),2)
    expect_equal(where_jk.cplist(a,3,5),3)
)

test_that("You can't find something that isn't in the cplist.", {
    expect_equal(where_jk.cplist(a,4,4),NA)
})



where_jk.cplist(a,4,4)



test_that("You aren't able to add (j,k,val) if (j,k) already exists in the matrix.", {
    expect_error(add.cplist(a,1,1,10/3))}
)


test_that("extract.cplist() works properly",{
    expect_equal(extract.cplist(a,1,1), 1.3)
})


test_that("extract.cplist() doesn't find something that's not in there",{
   expect_error(extract.cplist(a,1,3))
   expect_error(extract.cplist(a,10,1))
})



a=cplist(4)
a = add.cplist(a,0,1,1/3)
a = add.cplist(a,1,1,2/3)
a = add.cplist(a,2,1,5/3)


test_that("You can find valid existing indices the cplist.", {
    expect_equal()
add.cplist(a,2,2,3)
)

context("test plot")

## test .splitChangepoints

test_that("it makes a list with x and y components",{
  res <- .splitChangepoints(100, c(1:4), c(25,50,75))
  
  expect_true(length(res) == 4)
  
  nam <- sapply(res, names)
  expect_true(all(as.vector(nam) == c("x","y","x","y","x","y","x","y")))
  
  expect_true(all(res[[1]]$x == 1:25))
  expect_true(all(res[[2]]$x == 26:50))
  expect_true(all(res[[3]]$x == 51:75))
  expect_true(all(res[[4]]$x == 76:100))
  
  expect_true(all(res[[1]]$y == 1))
  expect_true(all(res[[2]]$y == 2))
  expect_true(all(res[[3]]$y == 3))
  expect_true(all(res[[4]]$y == 4))
  
  len <- sapply(res, function(x){length(x$y)})
  expect_true(all(len == 25))
})
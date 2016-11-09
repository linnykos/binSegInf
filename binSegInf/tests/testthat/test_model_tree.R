context("Test model tree")

## .create_node is correct

test_that(".create_node makes a valid node", {
  res <- .create_node(2, 5)
  expect_true(isValid(res))
})

test_that(".create_node errors if start is larger than end", {
  res <- expect_error(.create_node(5,2))
})

test_that(".create_node creates the right name", {
  res <- .create_node(2, 5)
  expect_true(res$name == "2-5")
})

##########################

## .get_leaves_names is correct

test_that(".get_leaves_names works", {
  data(acme, package = "data.tree")
  res <- .get_leaves_names(acme)
  expect_true(all(res == c("Go agile", "New Accounting Standards", "New Labs",
    "New Product Line", "New Software", "Outsource", "Switch to R")))
})

###############################

## .find_leadingBreakpoint is correct

test_that(".find_leadingBreakpoint will work if there is only one leaf", {
  tree <- .create_node(1, 10)
  tree$cusum <- 10
  
  res <- .find_leadingBreakpoint(tree)
  expect_true(res == "1-10")
})

#####################################

## .split_node is correct

test_that(".split_node correctly returns left and right", {
  tree <- .create_node(1, 10, 5)
  res <- .split_node(tree)
  
  expect_true(length(res) == 2)
  expect_true(res$left$start == 1)
  expect_true(res$left$end == 5)
  expect_true(res$right$start == 6)
  expect_true(res$right$end == 10)
})

#######################################

## .enumerate_splits is correct

test_that(".enumerate_splits is correct", {
  y <- c(rep(0, 10), rep(5, 5), rep(6, 5))
  res <- binSeg_fixedSteps(y, 2)
  
  expect_true(all(.enumerate_splits(res$tree) == c("1-20", "11-20")))
})
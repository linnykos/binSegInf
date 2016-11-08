context("Test model tree")

## .create_node is correct

test_that(".create_node makes a valid node", {
  res <- .create_node(2, 5)
  expect_true(isValid(res))
})

test_that(".create_node errors if start is larger than end", {
  res <- expect_error(.create_node(5,2))
})

##########################

## .get_leaves_names is correct

test_that(".get_leaves_names works", {
  data(acme, package = "data.tree")
  res <- .get_leaves_names(acme)
  expect_true(all(res == c("Go agile", "New Accounting Standards", "New Labs",
    "New Product Line", "New Software", "Outsource", "Switch to R")))
})
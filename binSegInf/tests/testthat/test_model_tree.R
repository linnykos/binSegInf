context("Test model tree")

## .create_node is correct

test_that(".create_node makes a valid node", {
  res <- .create_node(2, 5)
  expect_true(isValid(res))
})

test_that(".create_node errors if start is larger than end", {
  res <- expect_error(.create_node(5,2))
})
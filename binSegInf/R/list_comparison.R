.list_comparison.bsFs <- function(obj){
  is_valid(obj)

  numSteps <- obj$numSteps
  comp.lis <- vector("list", numSteps)
  active.vec <- .enumerate_splits(obj$tree)
  tree2 <- .create_node(1, obj$tree$end)

  for(i in 1:numSteps){
    breakpoint <- data.tree::FindNode(obj$tree, active.vec[i])$breakpoint
    comp.lis[[i]] <- .form_comparison(tree2, active.vec[i], breakpoint)

    node.selected <- data.tree::FindNode(tree2, active.vec[i])
    node.selected$breakpoint <- breakpoint
    
    node.pairs <- .split_node(node.selected)
    node.selected$AddChildNode(node.pairs$left)
    node.selected$AddChildNode(node.pairs$right)
  }

  comp.lis
}

.list_comparison.bsFs <- function(obj){
  is_valid(obj)

  numSteps <- obj$numSteps
  comp_lis <- vector("list", numSteps)
  active_vec <- .enumerate_splits(obj$tree)
  tree2 <- .create_node(1, obj$tree$end)

  for(i in 1:numSteps){
    breakpoint <- data.tree::FindNode(obj$tree, active_vec[i])$breakpoint
    comp_lis[[i]] <- .form_comparison(tree2, active_vec[i], breakpoint)

    node_selected <- data.tree::FindNode(tree2, active_vec[i])
    node_selected$breakpoint <- breakpoint
    
    node_pairs <- .split_node(node_selected)
    node_selected$AddChildNode(node_pairs$left)
    node_selected$AddChildNode(node_pairs$right)
  }

  comp_lis
}

.list_comparison.cbsFs <- function(obj){
  numSteps <- obj$numSteps
  comp_lis <- vector("list", numSteps)
  active_vec <- .enumerate_splits(obj$tree)
  tree2 <- .create_node(1, obj$tree$end)
  
  for(i in 1:numSteps){
    breakpoint <- data.tree::FindNode(obj$tree, active_vec[i])$breakpoint
    comp_lis[[i]] <- .form_comparison_cbs(tree2, active_vec[i], breakpoint)
    
    node_selected <- data.tree::FindNode(tree2, active_vec[i])
    node_selected$breakpoint <- breakpoint
    
    node_pairs <- .split_node_cbs(node_selected)
    if(!any(is.na(node_pairs$left))) node_selected$AddChildNode(node_pairs$left)
    node_selected$AddChildNode(node_pairs$middle)
    if(!any(is.na(node_pairs$right))) node_selected$AddChildNode(node_pairs$right)
  }
  
  comp_lis
}



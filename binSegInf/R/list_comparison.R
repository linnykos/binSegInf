.comparison_mat <- function(winning, losing){
  if(length(winning) != 3 | !is.numeric(winning)) stop("winning must be numeric of length 3")
  if(!is.matrix(losing) | !is.numeric(losing) | ncol(losing) != 3)
    stop("losing must be a numeric matrix with 3 columns")
  
  if(any(losing[,1] > losing[,2]) | any(losing[,2] >= losing[,3]))
    stop("losing must have second column >= first column and < third column")
  
  list(winning = winning, losing = losing)
}

# .list_comparison <- function(obj){
#   if(class(obj) != "bsFs") stop("obj must be class bsFs")
#   isValid(obj)
#   
#   numSteps <- obj$numSteps
#   comp.lis <- vector("list", numSteps)
#   active.vec <- .enumerate_splits(obj$tree)
#   tree2 <- .create_node(1, obj$tree$end)
#   
#   for(i in 1:numSteps){
#     comp.lis[[i]] <- .form_comparison(tree2, active.vec[i], 
#       obj$tree$FindNode(active.vec[i])$breakpoint)
#     
#     node.pairs <- .prepare_split(obj$tree$FindNode(active.vec[i]))
#     
#     node.selected <- tree2$FindNode(active.vec[i])
#     node.selected$AddChildNode(node.pairs$left)
#     node.selected$AddChildNode(node.pairs$right)
#   }
#   
#   comp.lis
# }
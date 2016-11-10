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

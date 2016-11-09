# binSeg_fixedSteps <- function(y, numSteps){
# 
#   #initialization
#   n <- length(y); tree <- .create_node(1, n)
# 
#   for(steps in 1:numSteps){
#     leaves.names <- .get_leaves_names(tree)
#     for(i in 1:length(leaves)){
#       res <- .find_breakpoint(y, tree$FindNode(leaves.names[i])$start, 
#         tree$FindNode(leaves.names[i])$end)
# 
#       tree$FindNode(leaves.names[i])$breakpoint <- res$breakpoint
#       tree$FindNode(leaves.names[i])$cusum <- res$cusum
#     }
#     
#     node.name <- .find_leadingBreakpoint(tree)
#   }
# 
#   structure(list(tree = tree, numSteps = numSteps))
# }

.find_leadingBreakpoint <- function(tree){
  leaves.names <- .get_leaves_names(tree)
  cusum.vec <- sapply(leaves.names, function(x){
    tree$FindNode(x)$cusum
  })
  leaves.names[which.max(cusum.vec)]
}

.find_breakpoint <- function(y, start, end){
  if(start >= end) stop("start must be smaller than end")
  
  idx <- seq(from = start, to = end - 1, by = 1)
  cusum.vec <- sapply(idx, .cusum, y = y, start = start, end = end)
  
  list(breakpoint = idx[which.max(cusum.vec)], cusum = max(cusum.vec))
}

.cusum <- function(y, start, idx, end){
  v <- .cusum_contrast(start, idx, end)
  v %*% y
}

#n1 is denoted as start to idx (inclusive)
.cusum_contrast <- function(start, idx, end){
  if(start > idx) stop("start must be smaller or equal than idx")
  if(idx >= end) stop("idx must be smaller to end")
  
  n1 <- idx - start + 1
  n2 <- end - idx
  
  c(rep(-1/n1, n1), rep(1/n2, n2)) * sqrt(1/((1/n1) + (1/n2)))
}
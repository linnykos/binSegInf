# binSeg_fixedSteps <- function(y, numSteps){
#   
#   #initialization
#   n <- length(y); tree <- .create_node(1, n)
#   
#   for(steps in 1:numSteps){
#     leaves.names <- .get_leaves_names(tree)
#     for(i in 1:length(leaves)){
#       res <- .find_breakpoint(y, population$FindNode(i)$start, 
#         population$FindNode(i)$end)
#       
#       population$FindNode(i)$breakpoint <- res$breakpoint
#       population$FindNode(i)$cusum <- res$cusum
#     }
#   }
#   
#   structure(list(tree = tree, numSteps = numSteps))
# }
# 
# .find_breakpoint <- function(y, start, end){
#   
# }
# 
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
  
  c(rep(-1/n1, n1), rep(1/n2, n2))
}
binSeg_fixedSteps <- function(y, numSteps){
    
  #initialization
  n <- length(y); tree <- .create_node(1, n)
    
  for(steps in 1:numSteps){
    leaves.names <- .get_leaves_names(tree)
    for(i in 1:length(leaves.names)){
      leaf <- tree$FindNode(leaves.names[i])
      
      res <- .find_breakpoint(y, leaf$start, leaf$end)
      
      leaf$breakpoint <- res$breakpoint; leaf$cusum <- res$cusum
    }
    
    node.name <- .find_leadingBreakpoint(tree)
    node.selected <- tree$FindNode(node.name)
    node.selected$active <- steps
     
    node.pairs <- .split_node(node.selected)
    node.selected$AddChildNode(node.pairs$left)
    node.selected$AddChildNode(node.pairs$right)
   }
    
  structure(list(tree = tree, numSteps = numSteps), class = "bsFs")
}

isValid.bsFs <- function(obj){
  if(class(obj$tree)[1] != "Node") stop("obj$tree must a Node")
  if(!is.numeric(obj$numSteps)) stop("obj$numSteps must be a numeric")
  if(length(.enumerate_splits(obj$tree)) != obj$numSteps) 
    stop("obj$tree and obj$numSteps disagree")
  
  TRUE
}

get_jumps.bsFs <- function(obj){
  get_jumps(obj$tree)
}

.find_breakpoint <- function(y, start, end){
  if(start > end) stop("start must be smaller than or equal to end")
  if(start == end) return(list(breakpoint = start, cusum = 0))
  
  idx <- seq(from = start, to = end - 1, by = 1)
  cusum.vec <- sapply(idx, .cusum, y = y, start = start, end = end)
  
  list(breakpoint = idx[which.max(abs(cusum.vec))], cusum = max(abs(cusum.vec)))
}

.cusum <- function(y, start, idx, end){
  v <- .cusum_contrast(start, idx, end)
  v %*% y[start:end]
}

#n1 is denoted as start to idx (inclusive)
.cusum_contrast <- function(start, idx, end){
  if(start > idx) stop("start must be smaller or equal than idx")
  if(idx >= end) stop("idx must be smaller to end")
  
  n1 <- idx - start + 1
  n2 <- end - idx
  
  c(rep(-1/n1, n1), rep(1/n2, n2)) * sqrt(1/((1/n1) + (1/n2)))
}

.cusum_contrast_full <- function(start, idx, end, n){
  res <- rep(0, n)
  res[start:end] <- .cusum_contrast(start, idx, end)
  
  res
}

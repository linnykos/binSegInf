circularBinSeg_fixedSteps <- function(y, numSteps){
  if(numSteps >= length(y)) stop("numSteps must be strictly smaller than the length of y")
  
  #initialization
  n <- length(y); tree <- .create_node(1, n); vec <- cumsum(y)

  for(steps in 1:numSteps){
    leaves.names <- .get_leaves_names(tree)
    for(i in 1:length(leaves.names)){
      leaf <- data.tree::FindNode(tree, leaves.names[i])
      
      res <- .find_breakpoint_cbs(vec[leaf$start:leaf$end])
      
      leaf$breakpoint <- res$breakpoint+leaf$start-1; leaf$cusum <- res$cusum
    }
    
    node.name <- .find_leadingBreakpoint(tree)
    node.selected <- data.tree::FindNode(tree, node.name)
    node.selected$active <- steps
    node.pairs <- .split_node_cbs(node.selected)
    if(!any(is.na(node.pairs$left))) node.selected$AddChildNode(node.pairs$left)
    node.selected$AddChildNode(node.pairs$middle)
    if(!any(is.na(node.pairs$right))) node.selected$AddChildNode(node.pairs$right)
  }
  
  obj <- structure(list(tree = tree, numSteps = numSteps), class = "cbsFs")
}

.find_breakpoint_cbs <- function(vec){
  n <- length(vec)
  
  breakpoint <- .enumerate_breakpoints_cbs(n)
  cusum_vec <- apply(breakpoint, 1, .cusum_cbs, vec = vec)
  
  idx <- which.max(abs(cusum_vec))
  list(breakpoint = breakpoint[idx,], cusum = cusum_vec[idx])
}

.enumerate_breakpoints_cbs <- function(n){
  cbind(rep(1:n, times = c((n-1), (n-1):1)), 
        unlist(lapply(1:n, function(x){
          if(x == 1) x:(n-1) else x:n
        })))
}

.cusum_cbs <- function(x, vec){
  stopifnot(length(x) == 2, x[1] >= 1, x[2] <= length(vec), x[1] <= x[2])
  stopifnot(all(x %% 1 == 0))
  stopifnot(!all(x[1]==1, x[2]==length(vec)))
  
  n <- length(vec); m <- x[2]-x[1]+1
  sum1 <- vec[x[2]] - ifelse(x[1] > 1, vec[x[1]-1], 0)
  sum2 <- (vec[n] - vec[x[2]]) + ifelse(x[1] > 1, vec[x[1]-1], 0)
  
  const <- sqrt(1/(1/m + 1/(n-m)))
  const*(sum1/m - sum2/(n-m))
}

.split_node_cbs <- function(node){
  if(any(is.na(node$breakpoint))) stop("node does not have a set breakpoint yet")
  stopifnot(length(node$breakpoint) == 2, node$breakpoint[1] >= node$start,
            node$breakpoint[2] <= node$end)
  
  if(node$breakpoint[1] > 1){left <- .create_node(node$start, node$breakpoint[1] - 1)} else{left <- NA}
  middle <- .create_node(node$breakpoint[1], node$breakpoint[2])
  if(node$breakpoint[2] < node$end){right <- .create_node(node$breakpoint[2] + 1, node$end)} else{right <- NA}
  
  stopifnot(!any(is.na(left)) | !any(is.na(right)))
  
  list(left = left, middle = middle, right = right)
}

is.na.Node <- function(x){FALSE}

.create_node <- function(start, end, breakpoint = NA, cusum = NA, active = NA){
  node <- data.tree::Node$new(paste0(start, "-", end))
  
  node$start <- start
  node$end <- end
  node$breakpoint <- breakpoint
  node$cusum <- cusum
  node$active <- active
  
  isValid(node)
  
  node
}

isValid.Node <- function(obj){
  if(obj$start > obj$end) stop("the start must be less or equal to end")
  if(!is.na(obj$breakpoint) & (obj$start > obj$breakpoint & obj$end < obj$breakpoint))
    stop("breakpoint must be between start and end (inclusive)")
  
  TRUE
}

.get_leaves_names <- function(tree){
  leaves <- tree$leaves
  vec <- sapply(leaves, function(x){x$name})
  names(vec) <- NULL
  sort(vec)
}

.find_leadingBreakpoint <- function(tree){
  leaves.names <- .get_leaves_names(tree)
  cusum.vec <- sapply(leaves.names, function(x){
    tree$FindNode(x)$cusum
  })
  leaves.names[which.max(abs(cusum.vec))]
}

.split_node <- function(node){
  if(is.na(node$breakpoint)) stop("node does not have a set breakpoint yet")
  if(node$breakpoint >= node$end) stop("node breakpoint must be less than end")
  
  left <- .create_node(node$start, node$breakpoint)
  right <- .create_node(node$breakpoint + 1, node$end)
  
  list(left = left, right = right)
}

.enumerate_splits <- function(tree){
  names(sort(tree$Get("active")))
}


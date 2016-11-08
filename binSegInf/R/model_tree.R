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

.get_leaves_names <- function(tree){
  leaves <- tree$leaves
  vec <- sapply(leaves, function(x){x$name})
  names(vec) <- NULL
  sort(vec)
}
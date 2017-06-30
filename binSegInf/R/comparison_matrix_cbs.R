.form_comparison_cbs <- function(tree, nodeName, breakpoint){
  stopifnot(length(breakpoint) == 2)
  
  winning <- matrix(NA, ncol = 4, nrow = 1)
  nodeNumeric <- .get_startEnd(nodeName)
  winning[c(1,4)] <- nodeNumeric; winning[2:3] <- breakpoint
  
  leaves.mat <- .get_leaves_matrix_excluding(tree, nodeName)
  
  losing <- .fourColumnMatrix_from_nodeVec(nodeNumeric, breakpoint)
  if(!any(is.na(leaves.mat))) {
    losing <- rbind(losing, do.call(rbind, .fourColumnMatrix_from_nodeMatrix(leaves.mat)))
  }
  
  list(winning = winning, losing = losing)
}

.fourColumnMatrix_from_nodeMatrix <- function(mat){
  plyr::alply(mat, 2, .fourColumnMatrix_from_nodeVec)
}

.fourColumnMatrix_from_nodeVec <- function(vec, exclude = NA){
  stopifnot(length(vec) == 2, is.numeric(vec))
  stopifnot(any(is.na(exclude)) || (is.numeric(exclude) & length(exclude) == 2) & 
              exclude[1] >= vec[1] & exclude[2] <= vec[2])
  
  if(any(is.na(exclude))){
    cbind(vec[1], .enumerate_breakpoints_cbs(vec[2]-vec[1]+1, vec[1]), vec[2])
  } else{
    mid.vec <- .enumerate_breakpoints_cbs(vec[2]-vec[1]+1, vec[1])
    mid.vec <- mid.vec[-which(apply(cbind(exclude[1] == mid.vec[,1], 
                                         exclude[2] == mid.vec[,2]), 1, all)),]
    if(length(mid.vec) == 0) return(numeric(0))
    cbind(vec[1], mid.vec, vec[2])
  }
}
polyhedra.cbsFt <- function(obj, ...){
  n <- .get_startEnd(obj$tree$name)[2] 
  comp_lis <- .list_comparison(obj)
  numNodes <- length(comp_lis)
  gamma_row_lis <- vector("list", numNodes)
  u_lis <- vector("list", numNodes)
  
  sign_vec <- .get_signs_cbsFt(obj)
  
  for(i in 1:numNodes){
    gamma_row_lis[[i]] <- .gammaRows_from_comparisons_cbsfs(comp_lis[[i]]$winning,
                                                            comp_lis[[i]]$losing, sign_vec[i], n)
    if(all(is.na(comp_lis[[i]]$winning))){
      stopifnot(nrow(gamma_row_lis[[i]]) %% 2 == 0)
      u_lis[[i]] <- c(rep(-obj$thres, nrow(gamma_row_lis[[i]])/2),
                      rep(-obj$thres, nrow(gamma_row_lis[[i]])/2))
    } else {
      u_lis[[i]] <- c(rep(0, nrow(gamma_row_lis[[i]]) - 1), obj$thres)
    }
  }
  
  mat <- do.call(rbind, gamma_row_lis)
  u <- do.call(c, u_lis)
  polyhedra(obj = mat, u = u)
}

.get_signs_cbsFt <- function(obj){
  nodes <- obj$tree$Get("active")
  active_vec <- names(sort(nodes))
  
  sign_vec <- rep(0, length(nodes))
  if(length(active_vec) == 0) return(sign_vec)
  
  for(i in 1:length(active_vec)){
    node_selected <- data.tree::FindNode(obj$tree, active_vec[i])
    sign_vec[i] <- sign(node_selected$cusum)
  }
  
  sign_vec
}
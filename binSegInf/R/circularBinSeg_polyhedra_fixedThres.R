polyhedra.cbsFt <- function(obj, ...){
  n <- .get_startEnd(obj$tree$name)[2] 
  comp_lis <- .list_comparison(obj)
  numJumps <- length(comp_lis)
  gamma_row_lis <- vector("list", numJumps)
  
  for(i in 1:numJumps){
    gamma_row_lis[[i]] <- .gammaRows_from_comparisons_cbsft(comp_lis[[i]]$winning,
                                                            comp_lis[[i]]$losing, sign_vec[i], n)
  }
  
  mat <- do.call(rbind, gamma_row_lis)
  u <- rep(0, nrow(mat))
  polyhedra(obj = mat, u = rep(0, nrow(mat)))
}



.gammaRows_from_comparisons_cbsft <- function(vec, mat, sign_win, n){
  
}
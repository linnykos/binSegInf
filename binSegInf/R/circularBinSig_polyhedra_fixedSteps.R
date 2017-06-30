polyhedra.cbsFs <- function(obj, ...){
  n <- .get_startEnd(obj$tree$name)[2] 
  numSteps <- obj$numSteps
  comp.lis <- .list_comparison(obj)
  sign.vec <- sign(jump_cusum(obj))
  gamma.row.lis <- vector("list", numSteps)
  
  for(i in 1:numSteps){
    losing.mat <- comp.lis[[i]]$losing
    
    gamma.row.lis[[i]] <- .gammaRows_from_comparisons_cbsfs(comp.lis[[i]]$winning,
                                                      losing.mat, sign.vec[i], n)
  }
  
  mat <- do.call(rbind, gamma.row.lis)
  polyhedra(obj = mat, u = rep(0, nrow(mat)))
}

.cusum_cbs_contrast_full <- function(start, idx, end, n){
  stopifnot(length(idx) == 2, idx[1] >= 1, idx[2] <= n, idx[1] <= idx[2])
  stopifnot(all(idx %% 1 == 0))
  stopifnot(!all(idx[1]==1, idx[2]==n))
  
  m <- idx[2]-idx[1]+1
  const <- sqrt(1/(1/m + 1/(n-m)))
  
  vec <- rep(0, n)
  if(idx[1] > start) vec[start:(idx[1]-1)] <- 1/(n-m) 
  vec[idx[1]:idx[2]] <- 1/m
  if(idx[2] < end) vec[(idx[2]+1):end] <- 1/(n-m)
  
  vec
}

.gammaRows_from_comparisons_cbsfs <- function(vec, mat, sign.win, n){
  stopifnot(length(vec) == 4, ncol(mat) == 4)
  
  win.contrast <- .cusum_cbs_contrast_full(vec[1], vec[2:3], vec[4], n)
  lose.contrast <- t(apply(mat, 1, function(x){
    .cusum_cbs_contrast_full(x[1], x[2:3], x[4], n)
  }))
  
  # add inequalities to compare winning split to all other splits
  res <- .vector_matrix_signedDiff(win.contrast, lose.contrast, sign.win, 
                                   rep(1, nrow(lose.contrast)))
  res2 <- .vector_matrix_signedDiff(win.contrast, lose.contrast, sign.win, 
                                    -rep(1, nrow(lose.contrast)))
  
  # add inequalities to compare splits to 0 (ensure correct sign)
  rbind(res, res2)
}
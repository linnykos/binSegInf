polyhedra.cbsFs <- function(obj, ...){
  n <- .get_startEnd(obj$tree$name)[2] 
  numSteps <- obj$numSteps
  comp_lis <- .list_comparison(obj)
  sign_vec <- sign(jump_cusum(obj))
  gamma_row_lis <- vector("list", numSteps)
  
  for(i in 1:numSteps){
    gamma_row_lis[[i]] <- .gammaRows_from_comparisons_cbsfs(comp_lis[[i]]$winning,
                                                            comp_lis[[i]]$losing, sign_vec[i], n)
  }
  
  mat <- do.call(rbind, gamma_row_lis)
  polyhedra(obj = mat, u = rep(0, nrow(mat)))
}

.cusum_cbs_contrast_full <- function(start, idx, end, n){
  stopifnot(length(idx) == 2, idx[1] >= 1, idx[2] <= n, idx[1] <= idx[2], start >= 1, end <= n)
  stopifnot(all(idx %% 1 == 0), start %% 1 == 0, end %% 1 == 0)
  stopifnot(!all(idx[1]==1, idx[2]==n))
  
  m <- idx[2]-idx[1]+1; n2 <- end-start+1
  const <- sqrt(1/(1/m + 1/(n2-m)))
  
  vec <- rep(0, n)
  if(idx[1] > start) vec[start:(idx[1]-1)] <- -1/(n2-m) 
  vec[idx[1]:idx[2]] <- 1/m
  if(idx[2] < end) vec[(idx[2]+1):end] <- -1/(n2-m)
  
  const * vec
}

.gammaRows_from_comparisons_cbsfs <- function(vec, mat, sign_win, n){
  stopifnot(length(vec) == 4, ncol(mat) == 4)
  
  win_contrast <- .cusum_cbs_contrast_full(vec[1], vec[2:3], vec[4], n)
  lose_contrast <- t(apply(mat, 1, function(x){
    .cusum_cbs_contrast_full(x[1], x[2:3], x[4], n)
  }))
  
  # add inequalities to compare winning split to all other splits
  res <- .vector_matrix_signedDiff(win_contrast, lose_contrast, sign_win, 
                                   rep(1, nrow(lose_contrast)))
  res2 <- .vector_matrix_signedDiff(win_contrast, lose_contrast, sign_win, 
                                    -rep(1, nrow(lose_contrast)))
  
  # add inequalities to compare splits to 0 (ensure correct sign)
  rbind(res, res2)
}
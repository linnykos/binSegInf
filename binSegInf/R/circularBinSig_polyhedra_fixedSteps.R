.cusum_cbs_contrast_full <- function(start, idx, end, n){
  stopifnot(length(idx) == 2, idx[1] >= 1, idx[2] <= n, idx[1] <= idx[2])
  stopifnot(all(idx %% 1 == 0))
  stopifnot(!all(idx[1]==1, idx[2]==length(vec)))
  
  m <- x[2]-x[1]+1
  const <- sqrt(1/(1/m + 1/(n-m)))
  
  vec <- rep(0, n)
  if(idx[1] > start) vec[start:(idx[1]-1)] <- 1/(n-m) 
  vec[idx[1]:idx[2]] <- 1/m
  if(idx[2] < end) vec[(idx[2]+1):end] <- 1/(n-m)
  
  vec
}
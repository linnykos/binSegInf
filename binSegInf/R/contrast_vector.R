contrast_vector.bsFs <- function(obj, jump.idx, sorted = F){
  jump.vec <- get_jumps(obj, sorted)
  jump <- jump.vec[jump.idx]
  
  n <- .get_startEnd(obj$tree$name)[2]
  jumpSorted.vec <- c(1, get_jumps(obj, T), n)
  idx <- which(jumpSorted.vec == jump)
  if(idx == 1) idx <- 2
  
  start <- jumpSorted.vec[idx-1]; split <- jump; end <- jumpSorted.vec[idx+1]
  res <- sign(.cusum_contrast_full(start, split, end, n))
  res[res>0] <- 1/sum(res > 0)
  res[res<0] <- -1/sum(res < 0)
  
  res
}
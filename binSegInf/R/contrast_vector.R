contrast_vector.bsFs <- function(obj, jump.idx, sorted = F, 
  type = c("segment", "spike")){
  if(!is.character(type[1])) stop("type must be character")
  if(!type[1] %in% c("segment", "spike")) stop("type must be either segment or spike")
  type <- type[1]
  
  jump.vec <- get_jumps(obj, sorted)
  jump <- jump.vec[jump.idx]
  n <- .get_startEnd(obj$tree$name)[2]
  
  if(type == "segment") {
    .contrast_vector_segment(obj, jump, n)
  } else {
    .contrast_vector_spike(obj, jump, n)
  }
}

.contrast_vector_segment <- function(obj, jump, n){
  jumpSorted.vec <- c(1, get_jumps(obj, T), n)
  
  idx <- which(jumpSorted.vec == jump)
  if(idx == 1) idx <- 2
  
  start <- jumpSorted.vec[idx-1]; split <- jump; end <- jumpSorted.vec[idx+1]
  res <- sign(.cusum_contrast_full(start, split, end, n))
  res[res>0] <- 1/sum(res > 0)
  res[res<0] <- -1/sum(res < 0)
  
  res
}

.contrast_vector_spike <- function(obj, jump, n){
  res <- rep(0, n)
  res[c(jump, jump+1)] <- c(-1, 1)
  
  res
}
.find_breakpoint_cbs <- function(vec){
  stopifnot(all(diff(vec) >= 0))
  
  n <- length(vec)
  
  breakpoint <- .enumerate_breakpoints_cbs(n)
  cusum_vec <- apply(breakpoint, 1, .cusum_cbs, vec = vec)
  
  idx <- which.max(abs(cusum_vec))
  list(breakpoint = breakpoint[idx,], cusum = cusum_vec[idx])
}

.enumerate_breakpoints_cbs <- function(n){
  cbind(rep(1:n, times = c((n-1), (n-1):1)), 
        unlist(lapply(1:n, function(x){
          if(x == 1) x:(n-1) else x:n
        })))
}

.cusum_cbs <- function(x, vec){
  stopifnot(all(diff(vec) >= 0))
  stopifnot(length(x) == 2, x[1] >= 1, x[2] <= length(vec), x[1] <= x[2])
  stopifnot(all(x %% 1 == 0))
  stopifnot(!all(x[1]==1, x[2]==length(vec)))
  
  n <- length(vec); m <- x[2]-x[1]+1
  sum1 <- vec[x[2]] - ifelse(x[1] > 1, vec[x[1]-1], 0)
  sum2 <- (vec[n] - vec[x[2]]) + ifelse(x[1] > 1, vec[x[1]-1], 0)
  
  const <- sqrt(1/m + 1/(n-m))
  const*(sum1/m - sum2/(n-m))
}

.split_node_cbs <- function(node){
  if(is.na(node$breakpoint)) stop("node does not have a set breakpoint yet")
  stopifnot(length(node$breakpoint) == 2, node$breakpoint[1] >= node$start,
            node$breakpoint[2] <= node$end)
  
  left <- ifelse(node$breakpoint[1] > 1, .create_node(node$start, node$breakpoint[1] - 1), NA)
  middle <- .create_node(node$breakpoint[1], node$breakpoint[2])
  right <- ifelse(node$breakpoint[2] < n, .create_node(node$breakpoint[2] + 1, node$end), NA)
  
  stopifnot(!is.na(left) | !is.na(right))
  
  list(left = left, middle = middle, right = right)
}
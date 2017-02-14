sample_splitting <- function(y, method, ...){
  n <- length(y)
  idx <- seq(2, n, by = 2)
  
  #estimate on half the data
  res <- jumps(method(y[idx], ...))
  
  #readjust the changepoints
  res <- (res-1)*2 + 1
  
  #average the points on the second half of the data
  structure(list(jumps = res, method = deparse(substitute(method))), class = "ss")
}

jumps.ss <- function(obj, sorted = F, ...){
  idx <- obj$jumps
  if(sorted) sort(idx) else idx
}

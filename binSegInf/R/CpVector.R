CpVector <- function(n, jump.height, jump.loc, func = stats::rnorm, ...){
  if(length(jump.height) != length(jump.loc) + 1) 
    stop("jump.height must be one element more than jump.loc")
  if(min(jump.loc) < 0 | max(jump.loc) >= 1)
    stop("jump.loc must lie between 0 (inclusive) and 1 (exclusive)")
  if(!all(diff(jump.loc) > 0)) stop("jump.loc must be strictly increasing")
  if(n %% 1 != 0 | n < 0) stop("n must be a positive integer")
  
  jump.idx <- .computeJumpIdx(n, jump.loc)
  mean.vec <- .formMeanVec(n, jump.height, jump.idx)
  
  data <- mean.vec + func(n, ...)
  
  structure(list(data = data, jump.height = jump.height, jump.idx = jump.idx),
    class = "CpVector")
}

.computeJumpIdx <- function(n, jump.loc){
  round(n*jump.loc)
}

.formMeanVec <- function(n, jump.height, jump.idx){
  diff.vec <- diff(c(0,jump.idx,n))
  rep(jump.height, times = diff.vec)
}

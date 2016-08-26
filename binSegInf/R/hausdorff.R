hausdorff <- function(set1, set2, one.sided = F){
  if(!is.numeric(set1) | !is.numeric(set2)) stop(paste("set1 and set2 must be",
    "numerics"))
  
  #handle corner cases
  if(length(set1) == 0 | length(set2) == 0) return(NA)
  if(length(set2) == 1) set2 = c(set2, set2)
  if(length(set1) == 1) set1 = c(set1, set1)

  dist.mat = sapply(set1, function(i){abs(i-set2)})
  if(class(dist.mat) != "matrix") dist.mat = as.matrix(dist.mat, nrow = 1)
  dist.vecx = apply(dist.mat, 2, min)
  
  if(!one.sided) dist.vecy = apply(dist.mat, 1, min) else dist.vecy = 0
  
  max(dist.vecx, dist.vecy)
}


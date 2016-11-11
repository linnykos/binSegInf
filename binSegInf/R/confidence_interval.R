confidence_interval <- function(y, polyhedra, contrast, sigma = 1, alpha = 0.05,
  gridsize = 250){
  seq.val <- seq(-2*min(abs(y)), 2*max(abs(y)), length.out = gridsize)
  pvalue <- pvalue(y, polyhedra, contrast, sigma, 
    value = seq.val)
  
  c(min(seq.val[pvalue >= alpha], na.rm = T), 
    max(seq.val[pvalue >= alpha], na.rm = T))
}
  
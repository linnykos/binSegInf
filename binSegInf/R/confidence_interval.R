#' Confidence interval from post-selection inference
#'
#' @param y numeric vector
#' @param polyhedra polyhedra object
#' @param contrast contrast numeric vector
#' @param sigma numeric to denote the sd of the residuals
#' @param alpha numeric between 0 and 1 with default of 0.05
#' @param gridsize numeric to denote how fine of a grid to invert the hypothesis
#' test
#'
#' @return a vector of two numbers, the lower and upper end of the confidence interval
#' @export
confidence_interval <- function(y, polyhedra, contrast, sigma = 1, alpha = 0.05,
  gridsize = 250){
  seq.val <- seq(-2*min(abs(y)), 2*max(abs(y)), length.out = gridsize)
  pvalue <- pvalue(y, polyhedra, contrast, sigma, 
    value = seq.val)
  
  c(min(seq.val[pvalue >= alpha], na.rm = T), 
    max(seq.val[pvalue >= alpha], na.rm = T))
}
  
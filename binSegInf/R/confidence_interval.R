#' Confidence interval from post-selection inference
#'
#' @param y numeric vector
#' @param polyhedra polyhedra object
#' @param contrast contrast numeric vector
#' @param sigma numeric to denote the sd of the residuals
#' @param alpha numeric between 0 and 1 with default of 0.05
#' @param gridsize numeric to denote how fine of a grid to invert the hypothesis
#' test
#' @param alternative string of either "one.sided" or "two.sided" for the 
#' alternative. If one.sided, the alternative means the test statistic is positive.
#'
#' @return a vector of two numbers, the lower and upper end of the confidence interval
#' @export
confidence_interval <- function(y, polyhedra, contrast, sigma = 1, alpha = 0.05,
  gridsize = 250, alternative = c("two.sided", "one.sided")){
  alternative <- match.arg(alternative, c("two.sided", "one.sided"))
  
  seq.val <- seq(-2*min(abs(y)), 2*max(abs(y)), length.out = gridsize)
  pvalue <- pvalue(y, polyhedra, contrast, sigma, null_mean = seq.val)
  
  if(alternative == "two.sided"){
    idx <- c(.select_index(pvalue, alpha, T), .select_index(pvalue, 1-alpha, F))
    c(seq.val[idx[1]], seq.val[idx[2]])
  } else {
    idx <- .select_index(pvalue, alpha, T)
    c(seq.val[idx], Inf)
  }
}
  
.select_index <- function(vec, alpha, lower = T){
  idx <- which.min(abs(vec - alpha))
  if(lower & vec[idx] > alpha & idx > 1) idx <- idx - 1
  if(!lower & vec[idx] < alpha & idx < length(vec)) idx <- idx + 1
  
  idx
}